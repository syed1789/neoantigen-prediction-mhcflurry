#!/usr/bin/env python3
"""
Fetch real Ensembl wildtype protein sequences for each unique protein_id
in a missense-variants JSON (produced by vep_annotate.py), and sanity
check that the reference amino acid VEP reported at each position
actually matches the fetched sequence -- if it doesn't, something is
wrong (offset bug, stale transcript version, etc.) and the variant is
dropped rather than silently mistranslated.
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.request

SEQ_URL = "https://rest.ensembl.org/sequence/id?type=protein"
BATCH_SIZE = 50
MAX_RETRIES = 5

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def to_one_letter(aa):
    """VEP's `amino_acids` field is already one-letter; guard against
    three-letter codes in case a caller feeds HGVSp text instead."""
    if len(aa) == 1:
        return aa
    return THREE_TO_ONE.get(aa, "X")


def fetch_batch(ids, timeout):
    payload = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request(
        SEQ_URL, data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read())


RETRYABLE_HTTP_CODES = {429, 500, 502, 503, 504}


def fetch_batch_with_retry(ids, timeout):
    delay = 2
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return fetch_batch(ids, timeout)
        except urllib.error.HTTPError as e:
            if e.code in RETRYABLE_HTTP_CODES and attempt < MAX_RETRIES:
                retry_after = float(e.headers.get("Retry-After", delay))
                time.sleep(retry_after)
                delay *= 2
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if attempt < MAX_RETRIES:
                time.sleep(delay)
                delay *= 2
                continue
            raise
    raise RuntimeError("unreachable")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("missense_json")
    ap.add_argument("--out", required=True, help="Output JSON: {protein_id: sequence}")
    ap.add_argument("--variants-out", required=True,
                     help="Filtered missense JSON, sanity-checked variants only")
    ap.add_argument("--sleep", type=float, default=0.4)
    ap.add_argument("--timeout", type=float, default=30.0)
    args = ap.parse_args()

    with open(args.missense_json) as f:
        variants = json.load(f)

    protein_ids = sorted({v["protein_id"] for v in variants if v.get("protein_id")})
    print(f"Fetching {len(protein_ids)} unique wildtype protein sequences from Ensembl",
          file=sys.stderr)

    sequences = {}
    for i in range(0, len(protein_ids), BATCH_SIZE):
        batch = protein_ids[i:i + BATCH_SIZE]
        print(f"  sequence batch {i // BATCH_SIZE + 1} "
              f"({i + 1}-{i + len(batch)} of {len(protein_ids)})", file=sys.stderr)
        try:
            results = fetch_batch_with_retry(batch, args.timeout)
        except Exception as e:
            print(f"  FAILED batch: {e} -- skipping", file=sys.stderr)
            continue
        for r in results:
            if "seq" in r:
                sequences[r["id"]] = r["seq"]
        time.sleep(args.sleep)

    print(f"Fetched {len(sequences)}/{len(protein_ids)} sequences", file=sys.stderr)

    kept, dropped = [], []
    for v in variants:
        pid = v.get("protein_id")
        pos = v.get("protein_position")
        ref_aa = to_one_letter(v.get("ref_aa", ""))
        seq = sequences.get(pid)
        if not seq or not pos:
            dropped.append((v, "no sequence or position"))
            continue
        if pos < 1 or pos > len(seq):
            dropped.append((v, "position out of range"))
            continue
        observed = seq[pos - 1]
        if observed != ref_aa:
            dropped.append((v, f"ref AA mismatch: VEP says {ref_aa}, "
                                f"fetched sequence has {observed} at position {pos}"))
            continue
        kept.append(v)

    print(f"Sanity check: {len(kept)} kept, {len(dropped)} dropped", file=sys.stderr)
    for v, reason in dropped:
        print(f"  dropped {v.get('gene_symbol')} {v.get('hgvsp')}: {reason}", file=sys.stderr)

    with open(args.out, "w") as f:
        json.dump(sequences, f)
    with open(args.variants_out, "w") as f:
        json.dump(kept, f, indent=2)


if __name__ == "__main__":
    main()
