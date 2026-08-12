#!/usr/bin/env nextflow

/*
 * Neoantigen candidate prediction, from a somatic VCF
 *
 * somatic VCF (PASS SNVs) -> Ensembl VEP REST annotation (missense only)
 *   -> real wildtype protein sequences (Ensembl REST)
 *   -> candidate 8-11mer peptides spanning each mutated residue
 *   -> MHCflurry MHC class I binding-affinity prediction
 *   -> ranked candidate neoantigen list
 *
 * Designed to consume the VCF output of ../nf-variant-pipeline's somatic.nf
 * (FilterMutectCalls' *.filtered.vcf.gz, or SnpEff's *.annotated.vcf.gz --
 * this pipeline does its own annotation and ignores any existing ANN field).
 *
 * Annotation uses the Ensembl VEP REST API rather than a local SnpEff/VEP
 * install: nf-variant-pipeline's own SNPEFF_ANNOTATE step documented that
 * the SnpEff database host was unreachable from its build environment (see
 * that repo's SOMATIC.md, "Known limitations"). VEP REST needs no database
 * download and was reachable here -- see NEOANTIGEN.md for what that
 * tradeoff costs (network dependency, rate limits) versus what it buys
 * (no multi-hundred-MB database, no host-availability roulette).
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------
params.vcf              = "$projectDir/test_data/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
params.vcf_filter        = "PASS"          // FILTER value to keep; the bundled HCC1395 truth VCF uses "PASS;HighConf"
params.chrom             = ""               // comma-separated chromosome restriction, e.g. "17,21"; empty = all
params.max_variants      = 0                // 0 = no cap; safety valve for ad-hoc runs against huge VCFs
params.hla_alleles       = "HLA-A*01:01,HLA-A*02:01,HLA-A*03:01,HLA-A*11:01,HLA-A*24:02,HLA-B*07:02,HLA-B*08:01,HLA-B*15:01,HLA-B*40:01,HLA-B*44:03,HLA-C*07:01,HLA-C*07:02"
params.peptide_lengths   = "8,9,10,11"
params.max_affinity_nm   = 500.0
params.top_n             = 50
// params.outdir default lives in nextflow.config

// sample_id used only for tagging/output filenames -- restricted to a safe
// charset for the same reason nf-variant-pipeline restricts sample_id
// (interpolated into filenames; an attacker-controlled VCF path shouldn't be
// able to inject shell metacharacters). See somatic.nf's validate_id().
def validate_id(id) {
    if (!(id ==~ /^[A-Za-z0-9_.-]+$/)) {
        error "Unsafe sample id '${id}' derived from --vcf filename — must match ^[A-Za-z0-9_.-]+\$"
    }
    id
}

// ---------------------------------------------------------------------
// ANNOTATE_VARIANTS — Ensembl VEP REST, canonical-transcript missense only
// ---------------------------------------------------------------------
process ANNOTATE_VARIANTS {
    tag "$sample_id"
    container 'python:3.11-slim@sha256:90744cff8f32887f075c47d747a173ff333e9e98801667af93c357fa9f5e28ff'
    publishDir "${params.outdir}/annotated", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    path script

    output:
    tuple val(sample_id), path("${sample_id}.missense_variants.json"), emit: variants

    script:
    def chrom_arg = params.chrom ? "--chrom ${params.chrom.split(',').join(' ')}" : ""
    def max_arg   = params.max_variants ? "--max-variants ${params.max_variants}" : ""
    """
    python3 "${script}" "${vcf}" \\
        --out "${sample_id}.missense_variants.json" \\
        --filter "${params.vcf_filter}" \\
        ${chrom_arg} ${max_arg}
    """
}

// ---------------------------------------------------------------------
// FETCH_PROTEIN_SEQUENCES — real Ensembl wildtype protein sequences,
// sanity-checked against VEP's reported reference amino acid
// ---------------------------------------------------------------------
process FETCH_PROTEIN_SEQUENCES {
    tag "$sample_id"
    container 'python:3.11-slim@sha256:90744cff8f32887f075c47d747a173ff333e9e98801667af93c357fa9f5e28ff'
    publishDir "${params.outdir}/annotated", mode: 'copy'

    input:
    tuple val(sample_id), path(variants_json)
    path script

    output:
    tuple val(sample_id), path("${sample_id}.wildtype_sequences.json"),
                           path("${sample_id}.verified_variants.json"), emit: sequences

    script:
    """
    python3 "${script}" "${variants_json}" \\
        --out "${sample_id}.wildtype_sequences.json" \\
        --variants-out "${sample_id}.verified_variants.json"
    """
}

// ---------------------------------------------------------------------
// GENERATE_PEPTIDES — 8-11mer windows spanning each mutated residue
// ---------------------------------------------------------------------
process GENERATE_PEPTIDES {
    tag "$sample_id"
    container 'python:3.11-slim@sha256:90744cff8f32887f075c47d747a173ff333e9e98801667af93c357fa9f5e28ff'
    publishDir "${params.outdir}/peptides", mode: 'copy'

    input:
    tuple val(sample_id), path(verified_variants_json), path(sequences_json)
    path script

    output:
    tuple val(sample_id), path("${sample_id}.candidate_peptides.tsv"), emit: peptides

    script:
    """
    python3 "${script}" "${verified_variants_json}" "${sequences_json}" \\
        --lengths "${params.peptide_lengths}" \\
        --out "${sample_id}.candidate_peptides.tsv"
    """
}

// ---------------------------------------------------------------------
// PREDICT_BINDING — MHCflurry class I binding-affinity prediction
// ---------------------------------------------------------------------
process PREDICT_BINDING {
    tag "$sample_id"
    container 'quay.io/biocontainers/mhcflurry:2.2.1--pyhdfd78af_0@sha256:3d6fe3090e69a904e7a77c11699d3948d8a76e3eb42e9a9b08ecf78647651dc2'
    publishDir "${params.outdir}/predictions", mode: 'copy'

    input:
    tuple val(sample_id), path(peptides_tsv)
    path script

    output:
    tuple val(sample_id), path("${sample_id}.raw_predictions.tsv"), emit: predictions

    script:
    """
    mhcflurry-downloads fetch models_class1_pan
    python3 "${script}" "${peptides_tsv}" \\
        --alleles "${params.hla_alleles}" \\
        --out "${sample_id}.raw_predictions.tsv"
    """
}

// ---------------------------------------------------------------------
// RANK_NEOANTIGENS — filter + sort by predicted affinity
// ---------------------------------------------------------------------
process RANK_NEOANTIGENS {
    tag "$sample_id"
    container 'quay.io/biocontainers/mhcflurry:2.2.1--pyhdfd78af_0@sha256:3d6fe3090e69a904e7a77c11699d3948d8a76e3eb42e9a9b08ecf78647651dc2'
    publishDir "${params.outdir}/predictions", mode: 'copy'

    input:
    tuple val(sample_id), path(predictions_tsv)
    path script

    output:
    path "${sample_id}.ranked_neoantigen_candidates.tsv"
    path "${sample_id}.top_candidates.md"

    script:
    """
    python3 "${script}" "${predictions_tsv}" \\
        --max-affinity-nm ${params.max_affinity_nm} \\
        --top-n ${params.top_n} \\
        --out "${sample_id}.ranked_neoantigen_candidates.tsv" \\
        --summary-out "${sample_id}.top_candidates.md"
    """
}

// ---------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------
workflow {
    log.info """
             NEOANTIGEN PREDICTION PIPELINE
             ===============================
             vcf              : ${params.vcf}
             vcf_filter       : ${params.vcf_filter}
             chrom            : ${params.chrom ?: '(all)'}
             hla_alleles      : ${params.hla_alleles}
             peptide_lengths  : ${params.peptide_lengths}
             max_affinity_nm  : ${params.max_affinity_nm}
             outdir           : ${params.outdir}
             """
             .stripIndent()

    sample_id = validate_id(file(params.vcf).simpleName)
    vcf_ch = Channel.of(tuple(sample_id, file(params.vcf, checkIfExists: true)))

    annotate_script = Channel.fromPath("${projectDir}/bin/vep_annotate.py").first()
    sequence_script  = Channel.fromPath("${projectDir}/bin/fetch_protein_sequences.py").first()
    peptide_script   = Channel.fromPath("${projectDir}/bin/generate_peptides.py").first()
    predict_script   = Channel.fromPath("${projectDir}/bin/run_mhcflurry_predict.py").first()
    rank_script      = Channel.fromPath("${projectDir}/bin/rank_neoantigens.py").first()

    ANNOTATE_VARIANTS(vcf_ch, annotate_script)
    FETCH_PROTEIN_SEQUENCES(ANNOTATE_VARIANTS.out.variants, sequence_script)

    peptide_input = FETCH_PROTEIN_SEQUENCES.out.sequences
        .map { id, seqs, verified -> tuple(id, verified, seqs) }
    GENERATE_PEPTIDES(peptide_input, peptide_script)

    PREDICT_BINDING(GENERATE_PEPTIDES.out.peptides, predict_script)
    RANK_NEOANTIGENS(PREDICT_BINDING.out.predictions, rank_script)
}
