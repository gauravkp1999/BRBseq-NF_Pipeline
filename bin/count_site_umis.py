#!/usr/bin/env python3
"""Count unique UMIs per polyA site and sample from STARsolo BAMs (reads carry CB/UB tags).

A read is assigned to the nearest polyA site lying 3' of its 3'-most aligned base (<= --window nt away, sense strand
only).  UMIs are de-duplicated per (sample, UMI, site).  Needs `samtools` on PATH.

  count_site_umis.py --pas reference/pas_reference.tsv --barcode-maps barcode_map_plate1.csv barcode_map_plate2.csv \
                     --bams plate1_Aligned.sortedByCoord.out.bam plate2_Aligned.sortedByCoord.out.bam -o site_counts.tsv
Per-sample BAMs (one sample each, no CB tag) can be given with --sample-bams name=path ...
"""
import argparse, bisect, csv, re, subprocess, sys
from collections import defaultdict

def load_pas(path):
    by = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by[(r["chrom"], r["strand"])].append((int(r["pos"]), r["pas_id"], r["gene_id"]))
    for k in by: by[k].sort()
    return by, {k: [x[0] for x in v] for k, v in by.items()}

def span(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cigar))

def stream(bam):
    p = subprocess.Popen(["samtools", "view", "-F", "4", bam], stdout=subprocess.PIPE, text=True)
    yield from p.stdout
    p.wait()

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pas", required=True)
    ap.add_argument("--bams", nargs="*", default=[])
    ap.add_argument("--barcode-maps", nargs="*", default=[], help="one per --bams entry, same order")
    ap.add_argument("--sample-bams", nargs="*", default=[], help="name=path for already-demultiplexed BAMs")
    ap.add_argument("--window", type=int, default=300)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    by, keys = load_pas(a.pas)
    umis = defaultdict(set)            # (site, sample) -> {umi}
    gene_of = {x[1]: x[2] for v in by.values() for x in v}
    jobs = []
    for bam, bm in zip(a.bams, a.barcode_maps):
        with open(bm) as f:
            m = {r["barcode"]: r["sample"] for r in csv.DictReader(f)}
        jobs.append((bam, m, None))
    for sb in a.sample_bams:
        n, p = sb.split("=", 1); jobs.append((p, None, n))
    samples = []
    for bam, m, fixed in jobs:
        for line in stream(bam):
            t = line.rstrip("\n").split("\t")
            flag, chrom, pos, cig = int(t[1]), t[2], int(t[3]), t[5]
            strand = "-" if flag & 16 else "+"
            tags = {x[:2]: x[5:] for x in t[11:]}
            ub = tags.get("UB", "-")
            if fixed: s = fixed
            else:
                s = m.get(tags.get("CB", "-"))
            if s is None or ub in ("-", ""): continue
            key = (chrom, strand)
            if key not in keys: continue
            arr = keys[key]
            if strand == "+":
                end = pos + span(cig) - 1
                i = bisect.bisect_left(arr, end - 3)
                if i < len(arr) and arr[i] - end <= a.window: site = by[key][i]
                else: continue
            else:
                i = bisect.bisect_right(arr, pos + 3) - 1
                if i >= 0 and pos - arr[i] <= a.window: site = by[key][i]
                else: continue
            umis[(site[1], s)].add(ub)
            if s not in samples: samples.append(s)
    samples = sorted(samples, key=lambda x: (x.split("_")[0] != "ctrl", x))
    with open(a.out, "w") as f:
        f.write("site_id\tgene_id\t" + "\t".join(samples) + "\n")
        for sid in sorted(gene_of):
            f.write(f"{sid}\t{gene_of[sid]}\t" + "\t".join(str(len(umis.get((sid, s), ()))) for s in samples) + "\n")

if __name__ == "__main__":
    main()
