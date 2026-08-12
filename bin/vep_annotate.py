#!/usr/bin/env python3
"""
Annotate a somatic VCF's SNVs against Ensembl VEP (REST) and emit the
canonical-transcript missense_variant hits as JSON.

Why REST and not a local SnpEff/VEP install: nf-variant-pipeline's own
SNPEFF_ANNOTATE step documented that the SnpEff database host was
unreachable from the build environment (see SOMATIC.md, "Known
limitations"). The Ensembl VEP REST API (https://rest.ensembl.org) needs
no database download and was reachable from this environment, so it is
used here as the annotation step instead. It is a real, versioned,
public annotation service (Ensembl VEP), not a local re-implementation.
"""
import argparse
import gzip
import json
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

VEP_URL = "https://rest.ensembl.org/vep/human/region?hgvs=1&canonical=1&protein=1"
BATCH_SIZE = 200
MAX_RETRIES = 5
# Ensembl's public REST server documents a ~15 req/s per-client budget;
# this stays well under that even with a full worker pool in flight.
DEFAULT_WORKERS = 6


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_snvs(vcf_path, filter_value, chrom_filter, max_variants):
    # Accept --chrom values with or without a "chr" prefix (e.g. "17" or
    # "chr17") regardless of whether the VCF itself uses one -- normalize
    # both sides before comparing. (Previously compared the raw VCF chrom,
    # e.g. "chr17", against a filter list of bare "17" and never matched
    # anything -- every chromosome-restricted run silently annotated zero
    # variants. See NEOANTIGEN.md debugging notes.)
    chrom_filter_norm = {c.replace("chr", "") for c in chrom_filter} if chrom_filter else None

    variants = []
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, _qual, filt = fields[:7]
            chrom_norm = chrom.replace("chr", "")
            if filter_value and filt != filter_value:
                continue
            if chrom_filter_norm and chrom_norm not in chrom_filter_norm:
                continue
            # SNVs only: this pipeline's peptide-generation step handles
            # single-residue substitutions, not indels/frameshifts.
            if len(ref) != 1 or len(alt) != 1 or ref == "-" or alt == "-":
                continue
            variants.append((chrom_norm, pos, ref, alt))
            if max_variants and len(variants) >= max_variants:
                break
    return variants


def vep_batch(batch, timeout):
    regions = [f"{c} {p} {p} {r}/{a} +" for c, p, r, a in batch]
    payload = json.dumps({"variants": regions}).encode()
    req = urllib.request.Request(
        VEP_URL,
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read())


RETRYABLE_HTTP_CODES = {429, 500, 502, 503, 504}


def vep_batch_with_retry(batch, timeout):
    delay = 2
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return vep_batch(batch, timeout)
        except urllib.error.HTTPError as e:
            if e.code in RETRYABLE_HTTP_CODES and attempt < MAX_RETRIES:
                retry_after = float(e.headers.get("Retry-After", delay))
                print(f"  VEP {e.code}, sleeping {retry_after}s "
                      f"(attempt {attempt}/{MAX_RETRIES})", file=sys.stderr)
                time.sleep(retry_after)
                delay *= 2
                continue
            raise
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < MAX_RETRIES:
                print(f"  VEP request error ({e}), retrying in {delay}s "
                      f"(attempt {attempt}/{MAX_RETRIES})", file=sys.stderr)
                time.sleep(delay)
                delay *= 2
                continue
            raise
    raise RuntimeError("unreachable")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf", help="Somatic VCF (optionally .gz), e.g. FilterMutectCalls output")
    ap.add_argument("--out", required=True, help="Output JSON path")
    ap.add_argument("--filter", default="PASS", help="Only keep records with this FILTER value "
                     "(default: PASS). Set to '' to disable filtering.")
    ap.add_argument("--chrom", nargs="*", default=None,
                     help="Restrict to these chromosomes (e.g. 17 chr17), default: all")
    ap.add_argument("--max-variants", type=int, default=0,
                     help="Cap number of SNVs read from the VCF, 0 = no cap")
    ap.add_argument("--sleep", type=float, default=1.0,
                     help="Seconds to sleep between VEP batches (politeness / rate limit; "
                     "applied per worker thread, not globally, when --workers > 1)")
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                     help="Concurrent VEP REST requests in flight (default: %(default)s)")
    args = ap.parse_args()

    variants = read_snvs(args.vcf, args.filter, args.chrom, args.max_variants)
    print(f"Read {len(variants)} candidate SNVs from {args.vcf} "
          f"(filter={args.filter!r}, chrom={args.chrom})", file=sys.stderr)

    batches = [variants[i:i + BATCH_SIZE] for i in range(0, len(variants), BATCH_SIZE)]
    n_batches = len(batches)
    missense = []
    missense_lock = threading.Lock()
    progress_lock = threading.Lock()
    done_count = [0]

    def process_batch(batch_num, batch):
        try:
            results = vep_batch_with_retry(batch, args.timeout)
        except Exception as e:
            print(f"  FAILED batch {batch_num}/{n_batches}: {e} -- skipping this batch",
                  file=sys.stderr)
            return
        local_hits = []
        for r in results:
            chrom, pos, _pos2, ref_alt, _strand = r.get("input", "").split()
            ref, alt = ref_alt.split("/")
            for tc in r.get("transcript_consequences", []):
                if tc.get("canonical") != 1:
                    continue
                if "missense_variant" not in tc.get("consequence_terms", []):
                    continue
                aa = tc.get("amino_acids", "")
                if "/" not in aa:
                    continue
                ref_aa, alt_aa = aa.split("/")
                local_hits.append({
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "gene_symbol": tc.get("gene_symbol"),
                    "gene_id": tc.get("gene_id"),
                    "transcript_id": tc.get("transcript_id"),
                    "protein_id": tc.get("protein_id"),
                    "hgvsp": tc.get("hgvsp"),
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "protein_position": tc.get("protein_start"),
                    "impact": tc.get("impact"),
                })
        with missense_lock:
            missense.extend(local_hits)
        with progress_lock:
            done_count[0] += 1
            if done_count[0] % 10 == 0 or done_count[0] == n_batches:
                print(f"  VEP batch {done_count[0]}/{n_batches} done "
                      f"({len(missense)} missense hits so far)", file=sys.stderr)
        if args.sleep:
            time.sleep(args.sleep)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(process_batch, i + 1, b) for i, b in enumerate(batches)]
        for f in as_completed(futures):
            f.result()

    print(f"Found {len(missense)} canonical-transcript missense variants", file=sys.stderr)
    with open(args.out, "w") as f:
        json.dump(missense, f, indent=2)


if __name__ == "__main__":
    main()
