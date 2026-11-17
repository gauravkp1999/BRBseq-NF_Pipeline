#!/usr/bin/env python3
"""Merge per-run gene x sample count matrices (TSV, gene_id as first column) into one combined gene x sample matrix."""
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Merge per-run BRB-seq count matrices into one matrix")
    parser.add_argument("matrices", nargs="+", help="Per-run gene x sample TSV files")
    parser.add_argument("-o", "--output", required=True, help="Output merged TSV file")
    args = parser.parse_args()

    merged = None
    seen_samples = set()
    for f in args.matrices:
        df = pd.read_csv(f, sep="\t", index_col=0)
        overlap = seen_samples & set(df.columns)
        if overlap:
            raise SystemExit(
                f"ERROR: sample name(s) {sorted(overlap)} appear in more than one run's "
                f"matrix. Sample names must be unique across all runs in the samplesheet."
            )
        seen_samples.update(df.columns)
        merged = df if merged is None else merged.join(df, how="outer")

    merged = merged.fillna(0).astype(int)
    merged.index.name = "gene_id"
    merged.to_csv(args.output, sep="\t")
    print(f"Merged {len(args.matrices)} run(s) into {args.output}: {merged.shape[0]} genes x {merged.shape[1]} samples")


if __name__ == "__main__":
    main()
