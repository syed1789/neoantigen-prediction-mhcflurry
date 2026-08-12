# Neoantigen prediction pipeline: methods, real results, and honest limitations

This extends the approach documented in
[nf-variant-pipeline](https://github.com/syed1789/nf-variant-pipeline)'s
`SOMATIC.md` to the next stage of a cancer genomics workflow: turning a
somatic call set into candidate MHC class I-binding peptides. Same
standard as that document: real data, real tool output, and an explicit
account of what broke and what the results do and don't mean.

## What this pipeline actually does

1. **Annotate**: send each somatic SNV to the [Ensembl VEP REST
   API](https://rest.ensembl.org/), keep only **canonical-transcript
   missense_variant** consequences (single amino acid substitutions).
2. **Fetch real protein sequences**: pull each affected transcript's
   wildtype protein sequence from Ensembl, and confirm the amino acid VEP
   reported at that position actually matches the fetched sequence before
   using it for anything (see "Debugging notes" below for why this check
   exists and what it catches).
3. **Generate peptides**: substitute the mutant residue into the wildtype
   sequence, then slice every 8/9/10/11-mer window that spans the mutated
   position.
4. **Predict binding**: run [MHCflurry](https://github.com/openvax/mhcflurry)
   2.2.1's `Class1AffinityPredictor` against a panel of common HLA class I
   alleles.
5. **Rank**: sort by predicted affinity (nM, ascending), keep predictions
   at or below 500nM (the standard "weak binder or better" convention),
   and compute a wildtype-vs-mutant fold change.

## What MHCflurry actually predicts — and what it doesn't

MHCflurry's `Class1AffinityPredictor` predicts **in vitro binding affinity
between a peptide and an MHC class I molecule**, in nanomolar IC50-scale
units, trained on quantitative binding assay and mass-spec data. That is
the *only* thing the number in `affinity_nm` means.

It does **not** predict, and this pipeline does **not** claim:

- **Antigen processing** — whether the proteasome would actually cut the
  source protein to produce this exact peptide (that's what MHCflurry's
  separate *processing* predictor / the combined *presentation* predictor
  models — not used here; see "Why binding-affinity-only" below).
- **Surface presentation** — whether the peptide is actually loaded onto
  MHC and displayed at the cell surface in this tumor.
- **TCR recognition or immunogenicity** — whether any T cell in the
  patient's repertoire would recognize this peptide-MHC complex at all.
  A peptide can bind MHC class I extremely well and be completely
  non-immunogenic (self-tolerance, no cognate TCR, etc.) — binding
  affinity is necessary but nowhere near sufficient for immunogenicity.
- **Expression** — whether the mutated transcript is even expressed in
  this tumor. This pipeline works purely from DNA-level variant calls; it
  has no RNA-seq input and makes no expression claim.
- **Confirmed genotype** — the panel of 12 HLA alleles used here (see
  below) is a common-population panel chosen for demonstration, **not**
  HCC1395's actual patient HLA genotype. That genotype was not reliably
  available to build this pipeline against (see "What I looked for and
  couldn't verify" below) — every allele referenced in this repo's example
  output should be read as "if a patient carried this allele," not as a
  claim about HCC1395's real MHC genotype.

A ranked row in `ranked_neoantigen_candidates.tsv` is a **computational
hypothesis about binding affinity only** — the same first-pass filter a
real neoantigen pipeline (pVACseq, and similar tools) applies before much
more expensive downstream validation (processing prediction, RNA
expression evidence, and ultimately wet-lab immunogenicity assays).
Presenting it as more than that would be exactly the kind of
scientifically misleading claim `SOMATIC.md` warned against when it
declined to run concordance benchmarking against mismatched data.

### Why binding-affinity-only, not the presentation predictor

MHCflurry also ships a `Class1PresentationPredictor` that combines binding
affinity with an antigen-processing model for a more realistic
presentation score. This pipeline uses the plain affinity predictor
because it's the well-established, most widely cited MHCflurry entry
point, and because scoping to one clearly-labeled prediction (affinity)
is more honest than blending two models into one number without being
able to validate the processing component's calibration in this
environment. Swapping in the presentation predictor is a natural, real
extension — noted here as future work, not silently done.

## What I looked for and couldn't verify: HCC1395's real HLA genotype

An early web search surfaced a specific-looking HLA-A/B/C genotype for
HCC1395. It did not check out: fetching the page it was attributed to
directly showed no HLA information on it at all, and no other source
checked corroborated it. Rather than use an unverified genotype and
present it as HCC1395's real MHC type, this pipeline uses a documented,
generic **common population-frequency HLA class I panel** instead (12
alleles: HLA-A\*01:01, A\*02:01, A\*03:01, A\*11:01, A\*24:02; B\*07:02,
B\*08:01, B\*15:01, B\*40:01, B\*44:03; C\*07:01, C\*07:02 — among the
most frequent class I alleles across global populations) and says so
explicitly, rather than fabricate a plausible-looking genotype. If
HCC1395's real HLA type is available from a source you trust, pass it via
`--hla_alleles` — the pipeline doesn't care where the allele list comes
from, it just shouldn't be presented as verified when it isn't.

## Real result: what genome-wide missense annotation of the HCC1395 truth set found

Full pipeline run against `test_data/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
(`--vcf_filter "PASS;HighConf"`), no chromosome restriction:

| Stage | Count |
|---|---|
| PASS;HighConf somatic SNVs in the truth VCF | 37,330 |
| Canonical-transcript missense_variant calls (VEP REST) | 294, across 284 distinct genes |
| Reference-AA sanity check (see debugging note 6) | 294/294 passed, 0 dropped |
| Unique wildtype protein sequences fetched (Ensembl REST) | 284/284 fetched |
| Candidate 8-11mer peptides generated | 10,957 |
| Peptide × allele affinity predictions (12-allele panel) | 131,484 |
| Predictions at or below 500nM ("weak binder or better") | 2,224 (209 strong ≤50nM, 2,015 weak ≤500nM), spanning 275 genes |

Committed as real, regeneratable output in `example_results/` — not hand-picked
or trimmed to look good; every row from that run is there.

Two concrete, real rows worth reading closely, for different reasons:

- **`TP53 p.Arg175His`** (`VRHCPHHER`, HLA-C\*07:01, 352.5nM, fold-change 4.1x vs.
  wildtype) — one of the most frequently reported TP53 mutations in cancer
  generally, found here as one independent hit among 294, not selected for or
  searched for in advance. A real missense call on a real, well-characterized
  cancer gene, correctly translated and correctly scored as a weak (not
  strong) predicted binder for this allele — MHCflurry does not treat "famous
  hotspot mutation" as a shortcut to a strong-binder prediction, which is the
  correct behavior: TP53's biological importance as a driver mutation and its
  MHC-I binding affinity are unrelated questions.
- **`CIC p.Pro895Leu`** (`GLASSPFQL`, HLA-A\*02:01, 24.7nM, fold-change ~190x
  vs. wildtype `GLASSPFQP` at 4,697.8nM) — the largest wildtype-to-mutant
  affinity swing near the top of the ranked list: the wildtype sequence is
  essentially a non-binder, and the single substitution turns it into a
  strong predicted binder. This is exactly the pattern `fold_change` is
  meant to surface — a mutation plausibly *creating* a new binding motif,
  as opposed to a mutation landing inside a peptide that would have bound
  well regardless (see rows with fold-change near 1.0x in the same table,
  e.g. `FIFTWCLQI`/MTRR at 0.999x).

## Debugging notes: real bugs found while building this

**1. MHCflurry 2.2.1 doesn't import on Python 3.13.**
`mhcflurry-downloads fetch` failed immediately with
`ModuleNotFoundError: No module named 'pipes'`. `pipes` was deprecated in
Python 3.11 and removed outright in 3.13; `mhcflurry/downloads_command.py`
in the latest published release (2.2.1, as of this writing) still imports
it unconditionally at module load time, so *any* use of the
`mhcflurry-downloads` CLI fails on 3.13, not just this pipeline's use of
it. Confirmed this isn't fixed by upgrading — 2.2.1 is the latest
non-release-candidate version on PyPI. Fixed by pinning the project to a
Python 3.11 virtualenv (`requirements.txt`, README) rather than patching
around a removed stdlib module — the correct fix is matching a supported
Python version, not carrying a shim indefinitely.

**2. `nf-variant-pipeline`'s SnpEff annotation limitation directly shaped
this pipeline's design, not just its documentation.** That repo's
`SOMATIC.md` documents that the SnpEff database host was unreachable from
its build environment. Rather than assume that was environment-specific
and re-attempt a local SnpEff/VEP install here, this pipeline was built
around the Ensembl VEP REST API from the start — which turned out to be
reachable. That's a real architecture decision driven by a real upstream
finding, not a hypothetical "you could also do X" note. The tradeoff is
explicit: REST annotation means this pipeline needs outbound network
access on every run and is subject to REST rate limits (handled here with
retry/backoff and a bounded thread pool — see `bin/vep_annotate.py`) —
where a correctly-cached local database wouldn't need either. For a
real institutional deployment processing patient data, a local,
periodically-refreshed VEP/SnpEff cache is very likely still the right
call; REST annotation of real patient variant coordinates against a
public third-party API is also a data-governance problem a clinical
pipeline can't wave away (see SECURITY.md).

**3. The initial sequential VEP annotation run was too slow to be
practical.** The first pass through all ~37,000 HCC1395 truth-set SNVs
sent one VEP batch request at a time with a 0.6s sleep between batches —
correct, but at an observed ~14 seconds per 200-variant batch (network +
server-side annotation time, not the sleep), a full genome-wide run
projected to roughly 45 minutes. Ensembl's REST server documents a much
higher per-client budget (order of 15 requests/second) than one
sequential request was using. Switched to a bounded thread pool (6
concurrent requests, `bin/vep_annotate.py --workers`) with the same
retry/backoff logic per request; the same genome-wide run (37,330 SNVs,
187 batches) completed in about 6 minutes instead of the ~45 minutes the
sequential version was on track for. Left the default worker count
conservative (6, not higher) rather than tuned to the exact edge of what
the server tolerates.

**4. The retry logic that made bug #3's concurrency safe had its own gap,
caught on the very next script.** `vep_batch_with_retry`'s exception
handler only retried on HTTP 429 (rate limiting); any other HTTP error
status re-raised immediately, no retry. This didn't surface during the
VEP annotation runs (no non-429 HTTP errors happened to occur), but the
identical pattern in `fetch_protein_sequences.py`'s retry helper did fail
immediately on the next run: fetching the 284 wildtype protein sequences
hit real, transient `500`/`503` responses from Ensembl on several
batches, each one immediately given up on rather than retried, and the
run finished "successfully" having fetched 34/284 sequences and dropped
all 294 variants for "no sequence." Fixed by widening the retryable-status
set to `{429, 500, 502, 503, 504}` in both scripts and re-running; the
corrected run fetched 284/284 sequences with no failures. Worth stating
plainly: this bug produced a **clean exit code and a real-looking (if
small) output file** — exactly the silent-wrong-answer failure mode this
document keeps returning to, not a crash that would have been hard to
miss.

**5. VEP REST silently omits `protein_id` unless you ask for it — caught
before it broke anything downstream, not after.** The first full
genome-wide annotation run completed cleanly (294 missense hits, no
errors) but every single hit had `protein_id: null`. VEP's `/vep/human/region`
endpoint only returns `protein_id` in `transcript_consequences` when the
request explicitly passes `protein=1` — `hgvs=1` alone gets you `hgvsp`
(a *string* like `ENSP00000269305.4:p.Arg175His`, useful for reading, not
for a REST call keyed on protein ID) but not the structured `protein_id`
field `fetch_protein_sequences.py` needs to actually fetch a sequence.
Caught by inspecting the raw JSON output before wiring it into the next
step, not by the pipeline failing loudly — feeding 294 null protein IDs
into `fetch_protein_sequences.py` would have fetched zero sequences and
silently produced an empty peptide list, the exact "reported success,
quietly did the wrong thing" failure mode `SOMATIC.md` warns about. Fixed
by adding `&protein=1` to `bin/vep_annotate.py`'s VEP REST query
(one-line fix; the genome-wide run was re-executed from scratch afterward
so every number in this document reflects the corrected query, not the
broken one).

**6. Reference-amino-acid mismatches are a real, not hypothetical, check —
even though this run didn't hit one.** `fetch_protein_sequences.py`
doesn't trust that a fetched Ensembl protein sequence at a VEP-reported
position is automatically the right one — it checks that the amino acid
VEP reported as the reference actually matches what the fetched sequence
has at that position, and drops (rather than silently mistranslates) any
variant where it doesn't. In this run, all 294/294 variants passed —
Ensembl's own VEP and sequence endpoints agreeing with each other isn't
surprising, since both ultimately read the same underlying annotation.
The check is still worth keeping and worth describing honestly as "ran
clean" rather than "caught a bug," because the two realistic failure
modes it guards against are real: a transcript whose *version* changed
between when VEP resolved the consequence and when the sequence endpoint
was queried, or (see debugging note 4) a batch that silently returned
partial/wrong data without raising an HTTP error at all. Silently
building a peptide around the wrong residue would look identical to
correct output unless something checks for it — the same class of
failure `SOMATIC.md` calls out for `nf-variant-pipeline`'s own silent
channel-exhaustion bugs.

**7. `--chrom` silently matched nothing, and the pipeline reported success
anyway.** `bin/vep_annotate.py --chrom` filters the VCF to specific
chromosomes before sending anything to VEP. The comparison used the raw
VCF chromosome string (`"chr17"`, since that's how this truth VCF names
its contigs) against the user-supplied filter list, which the README and
`neoantigen.nf` both document as bare numbers (`--chrom "17,21"`, no
`chr` prefix) — so `"chr17" not in ["17"]` was always true, and every
`--chrom`-restricted run silently annotated **zero** variants and
produced an empty (but well-formed) output all the way through the
pipeline, `[SUCCESS] completed=5 failed=0` and all. Caught by deliberately
testing `neoantigen.nf` end to end (not just the individual `bin/`
scripts) on a chromosome expected to contain a real hit
(chr17, which has the TP53 R175H call from the main run) and getting zero
results back instead. Fixed by normalizing the `chr` prefix off both the
VCF's chromosome value and the filter list before comparing, in
`read_snvs()`. Re-running the same chr17 restriction afterward recovered
exactly the same TP53 row (peptide, allele, affinity, fold-change) that
the unrestricted genome-wide run found — the fix didn't just stop
erroring, it reproduced the correct real answer.

This bug (and bugs #4 and #5 above) is why this document distinguishes
"ran without a crash" from "produced correct output" throughout — both
would have shipped with a passing-looking `[SUCCESS]` line if the
individual pipeline stages hadn't been re-run against a case whose
correct answer was already known in advance.

## Known limitations

- **Missense SNVs only.** This pipeline does not call or translate
  indels, frameshifts, splice-site variants, or gene fusions into
  peptides. Those variant classes require different translation logic
  (frameshift reading-frame tracking, splice-junction modeling) that
  isn't implemented here. `bin/vep_annotate.py` filters to single-base
  substitutions at the VCF-reading step specifically so this limitation
  fails loudly (variants silently dropped, visible in the read-vs-kept
  counts) rather than being fed into logic that assumes a single-residue
  substitution and produces a wrong answer.
- **Canonical transcript only.** A missense call on a non-canonical
  transcript is not surfaced, even if that transcript is the
  biologically relevant one in a given tissue. A production pipeline
  would need isoform-expression evidence (RNA-seq) to pick the right
  transcript per variant rather than always defaulting to canonical.
- **No expression filter.** Every missense variant that passes annotation
  becomes peptide candidates, regardless of whether the gene is expressed
  in this tumor. Real neoantigen pipelines cross-reference RNA-seq
  (variant-supporting read counts, or at minimum gene-level TPM) before
  prioritizing a candidate — omitted here because this pipeline's input
  is DNA-only, matching `nf-variant-pipeline`'s own scope.
- **Generic HLA panel, not a real patient genotype** — see above.
- **REST annotation dependency** — this pipeline requires outbound network
  access to `rest.ensembl.org` on every run; it is not a fully offline
  pipeline the way a pre-cached local SnpEff/VEP database would be.

## Talking points this gives you for an interview

- You can explain precisely what a binding-affinity number does and
  doesn't tell you about a "neoantigen candidate," and name the
  additional evidence (processing prediction, expression, immunogenicity
  assay) a real pipeline stacks on top before calling something
  validated.
- You can explain the wildtype-vs-mutant fold-change check and why
  ranking on raw mutant affinity alone is a recognized pitfall.
- You can speak to a real cross-repo architecture decision: choosing REST
  annotation here specifically *because* of a documented local-database
  failure in the sibling pipeline, including the tradeoffs that decision
  costs (network dependency, rate limits, data governance for real patient
  data).
- You can name one specific dependency/Python-version incompatibility you
  found, diagnosed, and fixed correctly (pinning the runtime, not patching
  around a removed module) — including reasoning about why it's the right
  fix relative to alternatives (shimming vs. version pinning).
