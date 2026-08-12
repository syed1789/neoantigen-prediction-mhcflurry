#!/usr/bin/env python3
"""
Turn sanity-checked missense variants + wildtype protein sequences into
candidate mutant peptides: every 8/9/10/11-mer window of the *mutant*
protein sequence that spans the substituted residue.

This is the standard "flanking window" approach used by neoantigen
pipelines for simple missense SNVs (e.g. pVACseq's equivalent step):
splice the single amino acid substitution into the wildtype sequence,
then slide each candidate MHC class I peptide length across the mutated
position. It does not attempt to model frameshifts, splice variants, or
fusions -- this pipeline only calls missense SNVs to begin with.
"""
import argparse
import csv
import json


def windows_spanning(seq, center_pos_1based, length):
    """All `length`-mers of `seq` that include position `center_pos_1based`
    (1-indexed). Returns 0-indexed (start, end_exclusive) tuples."""
    n = len(seq)
    center0 = center_pos_1based - 1
    out = []
    start_min = max(0, center0 - length + 1)
    start_max = min(center0, n - length)
    for start in range(start_min, start_max + 1):
        out.append((start, start + length))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("variants_json")
    ap.add_argument("sequences_json")
    ap.add_argument("--lengths", default="8,9,10,11")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.variants_json) as f:
        variants = json.load(f)
    with open(args.sequences_json) as f:
        sequences = json.load(f)
    lengths = [int(x) for x in args.lengths.split(",")]

    rows = []
    seen_peptides = set()
    for v in variants:
        wt_seq = sequences.get(v["protein_id"])
        if not wt_seq:
            continue
        pos = v["protein_position"]
        mut_seq = wt_seq[:pos - 1] + v["alt_aa"] + wt_seq[pos:]
        wt_peptide_full = wt_seq[max(0, pos - 12):pos + 11]  # context only, not used for prediction

        for length in lengths:
            for start, end in windows_spanning(mut_seq, pos, length):
                peptide = mut_seq[start:end]
                if "X" in peptide or "*" in peptide:
                    continue
                wt_window = wt_seq[start:end] if end <= len(wt_seq) else ""
                key = (peptide, v["hgvsp"])
                if key in seen_peptides:
                    continue
                seen_peptides.add(key)
                rows.append({
                    "peptide": peptide,
                    "length": length,
                    "gene_symbol": v["gene_symbol"],
                    "transcript_id": v["transcript_id"],
                    "protein_id": v["protein_id"],
                    "hgvsp": v["hgvsp"],
                    "chrom": v["chrom"],
                    "pos": v["pos"],
                    "ref": v["ref"],
                    "alt": v["alt"],
                    "protein_position": pos,
                    "ref_aa": v["ref_aa"],
                    "alt_aa": v["alt_aa"],
                    "wildtype_peptide": wt_window,
                })

    with open(args.out, "w", newline="") as f:
        fieldnames = list(rows[0].keys()) if rows else [
            "peptide", "length", "gene_symbol", "transcript_id", "protein_id",
            "hgvsp", "chrom", "pos", "ref", "alt", "protein_position",
            "ref_aa", "alt_aa", "wildtype_peptide",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} candidate peptides from {len(variants)} missense variants "
          f"(lengths {lengths})")


if __name__ == "__main__":
    main()
