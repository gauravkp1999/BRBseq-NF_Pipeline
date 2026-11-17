#!/usr/bin/env python3
"""
Convert a STARsolo raw/ matrix (matrix.mtx + barcodes.tsv + features.tsv) into a
gene x sample TSV count matrix, keeping only barcodes present in the barcode/sample
map (i.e. barcodes you actually used), and relabelling columns from raw barcode
sequences to human-readable sample names.

Parses the MatrixMarket (.mtx) coordinate format directly with pandas (no scipy
dependency) - STARsolo's matrix.mtx is a plain
    %%MatrixMarket matrix coordinate integer general
    % comment line(s)
    n_genes n_barcodes n_entries
    gene_idx barcode_idx count
    ...
(1-based indices, rows = genes/features, columns = barcodes).
"""
import argparse
import sys

import numpy as np
import pandas as pd


def read_mtx_triplets(path):
    """Return (n_rows, n_cols, DataFrame[row, col, value]) from a MatrixMarket file."""
    with open(path) as fh:
        header_seen = False
        dims = None
        for i, line in enumerate(fh):
            if line.startswith("%"):
                continue
            dims = tuple(int(x) for x in line.split())
            data_start = i + 1
            break
        else:
            sys.exit(f"ERROR: {path} has no data - could not find dimension line")

    n_rows, n_cols, _n_entries = dims
    triplets = pd.read_csv(
        path, sep=r"\s+", comment="%", skiprows=data_start, header=None,
        names=["row", "col", "value"], engine="python",
    )
    return n_rows, n_cols, triplets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, help="STARsolo matrix.mtx")
    parser.add_argument("--barcodes", required=True, help="STARsolo barcodes.tsv")
    parser.add_argument("--features", required=True, help="STARsolo features.tsv")
    parser.add_argument("--barcode-map", required=True, help="barcode<TAB>sample, no header")
    parser.add_argument("--run", required=True, help="Run/library ID, used in log messages")
    parser.add_argument("--min-reads", type=int, default=0, help="Drop samples with fewer total counts than this")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    barcodes = [line.strip() for line in open(args.barcodes) if line.strip()]
    features = pd.read_csv(args.features, sep="\t", header=None)
    gene_ids = features[0].astype(str).values

    # STARsolo matrix.mtx is features (rows) x barcodes (columns), 1-based indices
    n_rows, n_cols, triplets = read_mtx_triplets(args.matrix)
    if (n_rows, n_cols) != (len(gene_ids), len(barcodes)):
        sys.exit(
            f"ERROR: matrix dimensions {(n_rows, n_cols)} do not match "
            f"{len(gene_ids)} features x {len(barcodes)} barcodes"
        )

    # Dense gene x barcode array, filled from the sparse triplets (0 elsewhere)
    arr = np.zeros((n_rows, n_cols), dtype=np.int64)
    if not triplets.empty:
        arr[triplets["row"].to_numpy() - 1, triplets["col"].to_numpy() - 1] = triplets["value"].to_numpy()
    df = pd.DataFrame(arr, index=gene_ids, columns=barcodes)

    barcode_to_sample = {}
    with open(args.barcode_map) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            barcode_to_sample[parts[0].upper()] = parts[1]

    keep_cols = [bc for bc in df.columns if bc.upper() in barcode_to_sample]
    missing = set(barcode_to_sample) - {bc.upper() for bc in keep_cols}
    if missing:
        print(
            f"WARNING: {len(missing)} barcode(s) from the map were not found among "
            f"STARsolo's output barcodes for run '{args.run}': {sorted(missing)}",
            file=sys.stderr,
        )
    if not keep_cols:
        sys.exit(f"ERROR: none of the mapped barcodes for run '{args.run}' were found in the matrix")

    df = df[keep_cols]
    df.columns = [barcode_to_sample[bc.upper()] for bc in keep_cols]

    if args.min_reads > 0:
        totals = df.sum(axis=0)
        low = totals[totals < args.min_reads].index.tolist()
        if low:
            print(
                f"Dropping {len(low)} sample(s) below --min-reads {args.min_reads}: {low}",
                file=sys.stderr,
            )
        df = df.drop(columns=low)

    df.index.name = "gene_id"
    df.to_csv(args.output, sep="\t")
    print(f"Wrote {df.shape[0]} genes x {df.shape[1]} samples to {args.output}")


if __name__ == "__main__":
    main()
