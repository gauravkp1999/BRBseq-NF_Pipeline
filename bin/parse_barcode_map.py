#!/usr/bin/env python3
"""
Parse a user-supplied barcode/sample mapping file into two normalized outputs:
  1. whitelist.txt           - one barcode per line (for STARsolo --soloCBwhitelist)
  2. barcode_sample_map.tsv  - barcode<TAB>sample, no header (for matrix relabelling)

Accepts comma- or tab-delimited input, with a header row containing a sample-name
column (matched case-insensitively against: sample, samplename, sample_name, name)
and a barcode column (matched against: barcode, bc, b1, cb).
"""
import argparse
import csv
import sys

SAMPLE_ALIASES = {"sample", "samplename", "sample_name", "name"}
BARCODE_ALIASES = {"barcode", "bc", "b1", "cb"}


def sniff_delimiter(path):
    with open(path, newline="") as fh:
        first_line = fh.readline()
    if "\t" in first_line:
        return "\t"
    if "," in first_line:
        return ","
    raise SystemExit(f"Could not detect a delimiter (comma or tab) in {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("barcode_map", help="User-supplied barcode/sample mapping file")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    delim = sniff_delimiter(args.barcode_map)

    with open(args.barcode_map, newline="") as fh:
        reader = csv.reader(fh, delimiter=delim)
        rows = [row for row in reader if row and any(cell.strip() for cell in row)]

    if not rows:
        sys.exit(f"ERROR: {args.barcode_map} is empty")

    header = [c.strip().lower() for c in rows[0]]
    sample_idx = next((i for i, c in enumerate(header) if c in SAMPLE_ALIASES), None)
    barcode_idx = next((i for i, c in enumerate(header) if c in BARCODE_ALIASES), None)

    if sample_idx is None or barcode_idx is None:
        sys.exit(
            f"ERROR: could not find both a sample column (one of {sorted(SAMPLE_ALIASES)}) "
            f"and a barcode column (one of {sorted(BARCODE_ALIASES)}) in the header of "
            f"{args.barcode_map}. Found header: {rows[0]}"
        )

    data_rows = rows[1:]
    if not data_rows:
        sys.exit(f"ERROR: {args.barcode_map} has a header but no data rows")

    seen_barcodes = {}
    seen_samples = {}
    pairs = []
    for row in data_rows:
        if len(row) <= max(sample_idx, barcode_idx):
            continue
        sample = row[sample_idx].strip()
        barcode = row[barcode_idx].strip().upper()
        if not sample or not barcode:
            continue
        if barcode in seen_barcodes:
            sys.exit(f"ERROR: barcode '{barcode}' appears more than once in {args.barcode_map}")
        if sample in seen_samples:
            sys.exit(f"ERROR: sample name '{sample}' appears more than once in {args.barcode_map}")
        seen_barcodes[barcode] = sample
        seen_samples[sample] = barcode
        pairs.append((barcode, sample))

    if not pairs:
        sys.exit(f"ERROR: no valid barcode/sample rows parsed from {args.barcode_map}")

    with open(f"{args.outdir}/whitelist.txt", "w") as fh:
        for barcode, _ in pairs:
            fh.write(barcode + "\n")

    with open(f"{args.outdir}/barcode_sample_map.tsv", "w") as fh:
        for barcode, sample in pairs:
            fh.write(f"{barcode}\t{sample}\n")

    print(f"Parsed {len(pairs)} barcode/sample pairs from {args.barcode_map}")


if __name__ == "__main__":
    main()
