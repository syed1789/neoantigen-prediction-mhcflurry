# Security review

A deliberate pass over this repo's actual risk surface, in the same spirit
as `nf-variant-pipeline`'s inline security notes: specific to what this
pipeline does, not a generic checklist.

This is a security review, not a correctness review — several real
*correctness* bugs (an HTTP retry gap that silently dropped data, a
`--chrom` filter that silently matched nothing) were also found while
building this pipeline. They're documented in
[NEOANTIGEN.md](NEOANTIGEN.md#debugging-notes-real-bugs-found-while-building-this)
rather than here, since none of them are exploitable or represent a
trust-boundary problem — they're wrong-answer bugs, not attacker-facing
ones. Worth reading together: several of them share the same root cause
this document keeps calling out — code that exits cleanly while quietly
doing the wrong thing.

## Scope

This pipeline: (1) reads a local VCF file, (2) sends variant coordinates
and, later, peptide sequences to two public REST APIs (Ensembl, and
GitHub releases for MHCflurry model weights), (3) runs Python scripts
inside pinned containers or a local venv, (4) writes TSV/JSON/Markdown
output to disk. No credentials, no user accounts, no web-facing service.

## Findings

### 1. Real patient data would leave the machine over the network — flagged prominently, not just here

`bin/vep_annotate.py` and `bin/fetch_protein_sequences.py` send every
variant's chromosome/position/ref/alt, and every mutated transcript's ID,
to `rest.ensembl.org` (a third party) as part of normal operation — that's
the whole point of using VEP REST instead of a local database (see
NEOANTIGEN.md for why). For the real data bundled in this repo — GIAB/
SEQC2's public, consented HCC1395/HCC1395BL benchmark cell-line calls —
that's a non-issue; this data is specifically released for exactly this
kind of use. **It would not be a non-issue for a real patient VCF.**
Sending real patient variant coordinates to a public third-party API
without a data-processing agreement, encryption-in-transit review, and
governance sign-off is a real compliance problem (HIPAA in a US clinical
setting, equivalent regimes elsewhere), not a hypothetical one. This is
called out in both README.md and NEOANTIGEN.md, not buried here only,
specifically so it isn't missed by someone about to point `--vcf` at real
data. **Mitigation for real clinical use**: swap the annotation step for
a local, pre-cached VEP/SnpEff installation (no outbound variant data at
all) before running this against anything but public reference data.
Traffic to Ensembl is HTTPS, which limits exposure to passive network
eavesdropping but does nothing about the trust relationship with Ensembl
itself or data-retention questions.

### 2. Container pinning — fixed

Both containers referenced by `neoantigen.nf` (`python:3.11-slim`,
`quay.io/biocontainers/mhcflurry:2.2.1--pyhdfd78af_0`) are pinned by
**digest**, not just by tag, so a tag being force-pushed or a base image
being silently rebuilt upstream can't change what actually runs without
this repo's own Nextflow file changing. This goes one step further than
`nf-variant-pipeline`'s own containers (which are pinned by version tag,
not digest) specifically because this repo's build process made the
digests easy to capture at write time; re-verify and re-pin these
periodically, since a digest pin means updates require deliberate action,
not "it'll just get the latest patch."

### 3. Dependency pinning — fixed

`requirements.txt` pins `mhcflurry==2.2.1` exactly (required — this is the
version whose Python 3.13 incompatibility was diagnosed and worked around;
see NEOANTIGEN.md) and pins `pandas`, `biopython`, and `requests` to the
exact versions actually tested (not floating `>=` ranges), for the same
supply-chain-predictability reason `nf-variant-pipeline` pins specific
container tags rather than `:latest`.

### 4. Shell interpolation in Nextflow `script:` blocks

Every `neoantigen.nf` process interpolates path variables (`"${vcf}"`,
`"${peptides_tsv}"`, etc.) directly into double-quoted shell command
strings — the same pattern `nf-variant-pipeline`'s `somatic.nf` and
`main.nf` use. Double-quoting neutralizes whitespace and most shell
metacharacters in a filename, but not `$`, backticks, or `\` — a
maliciously crafted input filename containing those could still inject
into the command. `nf-variant-pipeline` addresses this specifically for
`sample_id` (validated with `^[A-Za-z0-9_.-]+$` before it's used to name
outputs or appear in read-group strings), since `sample_id` is a string
value under direct control in that pipeline. This pipeline follows the
same convention: `neoantigen.nf`'s `validate_id()` restricts the
Nextflow-derived `sample_id` to the same safe charset before it's used
in any output filename. The one thing neither pipeline sanitizes is the
*raw filename* of a `path`-typed input itself (`params.vcf`, `params.reads`
upstream) — Nextflow stages those files into the work directory under
their original basename, and that basename is what gets interpolated.
Exploiting this requires an attacker to control the literal filename
passed via `--vcf`, which is a much narrower bar than arbitrary shell
metacharacter injection via unquoted interpolation, and is accepted here
consistently with the upstream repo's existing risk posture rather than
introduced fresh. Documented rather than silently left unmentioned.

### 5. No secrets, no credentials, no persistent state

Nothing in this pipeline reads or writes credentials, API keys, or
tokens. The Ensembl and GitHub-releases endpoints used are unauthenticated
public APIs. No database, no user accounts, nothing to leak beyond the
input VCF's contents (see finding 1) and the pipeline's own output.

### 6. Input validation on the VCF itself

`bin/vep_annotate.py` parses the VCF with plain `str.split("\t")`, not a
VCF-format-aware library — sufficient for well-formed VCFs (real Mutect2 /
GATK / SomaticSeq output, the only realistic input for this pipeline) but
not hardened against a deliberately malformed file. This is a correctness
concern more than a security one (a malformed line just fails to unpack
and raises `ValueError`, it doesn't execute anything), but is noted since
"parses untrusted structured input with manual string splitting" is
usually worth a second look. No fix applied — a full VCF parser
(`pysam`/`cyvcf2`) would be the right call for a pipeline meant to run
untrusted third-party VCFs; this one is scoped to well-formed
Mutect2-family output, consistent with the rest of this repo's scope.

## Summary

| # | Finding | Severity | Status |
|---|---|---|---|
| 1 | Real patient data would transit a third-party API if misused | High (conditional — not applicable to the public data this repo ships with) | Documented prominently in README, NEOANTIGEN.md, and here; not fixable in-pipeline (it's an intended tradeoff for public data, wrong tradeoff for patient data) |
| 2 | Containers pinned by tag only | Low | Fixed — pinned by digest |
| 3 | Dependencies pinned to loose ranges | Low | Fixed — pinned to tested versions |
| 4 | Shell interpolation of path basenames | Low, same posture as upstream repo | Documented, consistent with `nf-variant-pipeline`'s accepted risk |
| 5 | Secrets handling | None found | N/A |
| 6 | Manual VCF parsing | Low (correctness, not security) | Documented, not fixed — out of current scope |
