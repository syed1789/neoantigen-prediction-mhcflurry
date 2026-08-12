#!/usr/bin/env python3
"""
Run MHCflurry (Class1AffinityPredictor, BA-only -- see NEOANTIGEN.md for
why not the presentation predictor) against a panel of HLA class I
alleles for every candidate peptide.

MHCflurry predicts *binding affinity* (nM, lower = tighter predicted
binding) between a peptide and an MHC class I allele. It does not
predict antigen processing, surface presentation, TCR recognition, or
immunogenicity -- see NEOANTIGEN.md, "What this pipeline does NOT tell
you."
"""
import argparse
import csv
import sys


DEFAULT_ALLELES = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*40:01", "HLA-B*44:03",
    "HLA-C*07:01", "HLA-C*07:02",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("peptides_tsv", help="Output of generate_peptides.py")
    ap.add_argument("--alleles", default=",".join(DEFAULT_ALLELES),
                     help="Comma-separated HLA class I alleles (MHCflurry allele-name format)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    alleles = [a.strip() for a in args.alleles.split(",") if a.strip()]

    with open(args.peptides_tsv, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    if not rows:
        print("No candidate peptides -- writing header-only output", file=sys.stderr)
        header = [
            "peptide", "length", "allele", "affinity_nm", "percentile_rank",
            "wildtype_peptide", "wildtype_affinity_nm",
            "gene_symbol", "hgvsp", "transcript_id", "chrom", "pos", "ref", "alt",
            "protein_position", "ref_aa", "alt_aa",
        ]
        with open(args.out, "w", newline="") as f:
            csv.writer(f, delimiter="\t").writerow(header)
        return

    # Import here so --help works without mhcflurry installed / models fetched.
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()

    supported = set(predictor.supported_alleles)
    unsupported = [a for a in alleles if a not in supported]
    if unsupported:
        print(f"WARNING: these alleles are not supported by the loaded MHCflurry "
              f"models and will be skipped: {unsupported}", file=sys.stderr)
    alleles = [a for a in alleles if a in supported]
    if not alleles:
        raise SystemExit("No requested alleles are supported by the loaded MHCflurry models.")

    peptides = [r["peptide"] for r in rows]
    peptide_lengths_ok = [8 <= len(p) <= 15 for p in peptides]
    if not all(peptide_lengths_ok):
        bad = sum(1 for ok in peptide_lengths_ok if not ok)
        print(f"WARNING: {bad} peptides outside MHCflurry's supported 8-15mer "
              f"range will be dropped by the predictor", file=sys.stderr)

    print(f"Predicting {len(peptides)} peptides x {len(alleles)} alleles "
          f"({len(peptides) * len(alleles)} peptide-allele pairs)", file=sys.stderr)

    all_predictions = []
    for allele in alleles:
        df = predictor.predict_to_dataframe(peptides=peptides, allele=allele, throw=False)
        df["allele"] = allele
        all_predictions.append(df)

    import pandas as pd
    predictions = pd.concat(all_predictions, ignore_index=True)

    meta = pd.DataFrame(rows)
    merged = predictions.merge(meta, on="peptide", how="left")

    merged = merged.rename(columns={
        "prediction": "affinity_nm",
        "prediction_percentile": "percentile_rank",
    })

    # Also predict the corresponding wildtype peptide's affinity for the same
    # allele, so downstream ranking can flag mutant-specific binders (a
    # mutant peptide that binds well only because the *wildtype* sequence
    # already binds well isn't evidence the mutation created anything new).
    wt_pairs = merged[["wildtype_peptide", "allele"]].dropna()
    wt_pairs = wt_pairs[wt_pairs["wildtype_peptide"].str.len().between(8, 15)]
    wt_pairs = wt_pairs.drop_duplicates()
    wt_by_allele = {}
    for allele, group in wt_pairs.groupby("allele"):
        wt_peptides = sorted(set(group["wildtype_peptide"]))
        wt_df = predictor.predict_to_dataframe(peptides=wt_peptides, allele=allele, throw=False)
        wt_by_allele[allele] = dict(zip(wt_df["peptide"], wt_df["prediction"]))

    def lookup_wt(row):
        return wt_by_allele.get(row["allele"], {}).get(row.get("wildtype_peptide"))

    merged["wildtype_affinity_nm"] = merged.apply(lookup_wt, axis=1)

    keep_cols = [
        "peptide", "length", "allele", "affinity_nm", "percentile_rank",
        "wildtype_peptide", "wildtype_affinity_nm",
        "gene_symbol", "hgvsp", "transcript_id", "chrom", "pos", "ref", "alt",
        "protein_position", "ref_aa", "alt_aa",
    ]
    keep_cols = [c for c in keep_cols if c in merged.columns]
    merged = merged[keep_cols]
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} predictions to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
