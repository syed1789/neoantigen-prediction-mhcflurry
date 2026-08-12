#!/usr/bin/env python3
"""
Collapse raw MHCflurry peptide x allele predictions into a ranked
candidate neoantigen list.

Ranking is by predicted mutant binding affinity (nM), ascending
(tighter predicted binding first) -- the same convention IEDB/NetMHC
tools use. Standard, widely-cited affinity bands (not a hard biological
cutoff -- see NEOANTIGEN.md) are added as a `binder_class` column:
strong <=50nM, weak <=500nM, non-binder above that.

`fold_change` = wildtype_affinity_nm / affinity_nm. A high fold_change
means the mutation substantially improved predicted binding relative to
the unmutated sequence -- the mutation, not a pre-existing motif, is
plausibly responsible for the predicted binder. A fold_change near 1
means the wildtype peptide was already predicted to bind about as well,
so the mutation isn't what's driving the (predicted) binding.
"""
import argparse
import sys

import pandas as pd


def binder_class(affinity_nm):
    if pd.isna(affinity_nm):
        return "NA"
    if affinity_nm <= 50:
        return "strong"
    if affinity_nm <= 500:
        return "weak"
    return "non-binder"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("predictions_tsv")
    ap.add_argument("--max-affinity-nm", type=float, default=500.0,
                     help="Only keep predictions at or below this affinity (default: 500nM, "
                     "the standard 'weak binder or better' threshold)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=50,
                     help="Also write this many top candidates as a Markdown summary table")
    ap.add_argument("--summary-out", default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.predictions_tsv, sep="\t")
    if df.empty:
        print("No predictions to rank -- writing empty ranked output "
              "and summary (both outputs are always written, even when empty, "
              "so a Nextflow run doesn't fail on a missing declared output).",
              file=sys.stderr)
        df.insert(0, "rank", [])
        df["binder_class"] = []
        df["fold_change_vs_wildtype"] = []
        df.to_csv(args.out, sep="\t", index=False)
        if args.summary_out:
            with open(args.summary_out, "w") as f:
                f.write("# Top candidate neoantigens (predicted MHC class I binding affinity)\n\n"
                        "No candidate peptides were generated for this input "
                        "(e.g. a chromosome/region restriction with no missense "
                        "hits) -- nothing to rank.\n")
        return

    df["binder_class"] = df["affinity_nm"].apply(binder_class)
    df["fold_change_vs_wildtype"] = df["wildtype_affinity_nm"] / df["affinity_nm"]

    ranked = df[df["affinity_nm"] <= args.max_affinity_nm].copy()
    ranked = ranked.sort_values("affinity_nm", ascending=True)
    ranked.insert(0, "rank", range(1, len(ranked) + 1))

    ranked.to_csv(args.out, sep="\t", index=False)
    print(f"{len(ranked)}/{len(df)} peptide-allele predictions at or below "
          f"{args.max_affinity_nm}nM written to {args.out}", file=sys.stderr)

    if args.summary_out:
        top = ranked.head(args.top_n)
        lines = [
            "# Top candidate neoantigens (predicted MHC class I binding affinity)",
            "",
            f"{len(ranked)} peptide/allele pairs at or below {args.max_affinity_nm}nM, "
            f"out of {len(df)} peptide/allele pairs scored.",
            "",
            "Ranked by predicted binding affinity only -- this is *not* a ranking of "
            "confirmed neoantigens. See NEOANTIGEN.md before treating any row here as "
            "more than a computational hypothesis.",
            "",
            "| Rank | Peptide | Allele | Affinity (nM) | Class | Gene | HGVSp | "
            "Fold change vs WT |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for _, r in top.iterrows():
            fc = r["fold_change_vs_wildtype"]
            fc_str = f"{fc:.1f}x" if pd.notna(fc) else "NA"
            lines.append(
                f"| {r['rank']} | {r['peptide']} | {r['allele']} | "
                f"{r['affinity_nm']:.1f} | {r['binder_class']} | {r['gene_symbol']} | "
                f"{r['hgvsp']} | {fc_str} |"
            )
        with open(args.summary_out, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"Top-{args.top_n} summary written to {args.summary_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
