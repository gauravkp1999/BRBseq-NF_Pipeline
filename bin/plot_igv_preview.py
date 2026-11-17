#!/usr/bin/env python3
"""Static 'IGV-like' preview: mean coverage (per million aligned reads) of ctrl vs KD over one gene per event type.
Genes are drawn 5'->3' (minus-strand genes flipped).  Needs samtools + matplotlib.
  plot_igv_preview.py DATA_DIR BAM_DIR out.png
"""
import csv, subprocess, sys, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D, B, out = sys.argv[1:4]
meta = {r["sample"]: r["condition"] for r in csv.DictReader(open(f"{D}/sample_meta.csv"))}
groups = {c: [f"{B}/{s}.bam" for s, k in meta.items() if k == c] for c in ("ctrl", "KD")}
tot = {b: int(subprocess.run(["samtools", "view", "-c", "-F", "4", b], capture_output=True, text=True).stdout) for g in groups.values() for b in g}
gtf = collections.defaultdict(list)
for l in open(f"{D}/reference/genes.gtf"):
    t = l.rstrip("\n").split("\t")
    gid = t[8].split('gene_id "')[1].split('"')[0]
    gtf[gid].append((t[2], t[0], int(t[3]), int(t[4]), t[6]))
pas = collections.defaultdict(list)
for r in csv.DictReader(open(f"{D}/reference/pas_reference.tsv"), delimiter="\t"): pas[r["gene_id"]].append(r)
tg = list(csv.DictReader(open(f"{D}/expected/truth_genes.tsv"), delimiter="\t"))
truth_sites = {r["site_id"]: r for r in csv.DictReader(open(f"{D}/expected/truth_sites.tsv"), delimiter="\t")}
pick = []
for cat in ("ipa_up", "ipa_down", "expr_only", "trap"):
    for r in tg:
        if r["event_type"] == cat and (cat != "trap" or float(truth_sites[r["tested_site_id"]]["delta_frac"]) > 0.1):
            pick.append(r); break
fig, ax = plt.subplots(3, len(pick), figsize=(4.3 * len(pick), 6.4), gridspec_kw={"height_ratios": [1, 3, 3]}, squeeze=False)
col = {"ctrl": "#3b6ea5", "KD": "#d9822b"}
for j, r in enumerate(pick):
    g = r["gene_id"]; feats = gtf[g]; gene = [f for f in feats if f[0] == "gene"][0]
    chrom, s, e, strand = gene[1], gene[2], gene[3], gene[4]
    x = np.arange(s, e + 1); flip = strand == "-"
    pos = (e - x) if flip else (x - s)
    for cond, bams in groups.items():
        cov = np.zeros(len(x))
        for b in bams:
            d = subprocess.run(["samtools", "depth", "-a", "-Q", "0", "-r", f"{chrom}:{s}-{e}", b], capture_output=True, text=True).stdout
            cov += np.array([int(l.split("\t")[2]) for l in d.splitlines()]) / tot[b] * 1e6
        cov /= len(bams)
        a = ax[1 if cond == "ctrl" else 2, j]
        a.fill_between(pos, cov, step="mid", color=col[cond], lw=0); a.set_ylabel(f"{cond}\nmean CPM" if j == 0 else "")
        a.set_xlim(0, len(x)); a.spines[["top", "right"]].set_visible(False)
        ax[1, j].set_ylim(0, None)
    top = max(ax[1, j].get_ylim()[1], ax[2, j].get_ylim()[1]); ax[1, j].set_ylim(0, top); ax[2, j].set_ylim(0, top)
    m = ax[0, j]; m.set_xlim(0, len(x)); m.set_ylim(0, 1); m.axis("off")
    m.plot([0, len(x)], [.5, .5], color="#888", lw=1)
    for f in feats:
        if f[0] == "exon":
            a0, a1 = sorted(((e - f[3]), (e - f[2])) if flip else ((f[2] - s), (f[3] - s)))
            m.add_patch(plt.Rectangle((a0, .25), a1 - a0, .5, color="#444"))
    for p in pas[g]:
        px = (e - int(p["pos"])) if flip else (int(p["pos"]) - s)
        for k in (1, 2): ax[k, j].axvline(px, color="#c0392b", lw=.8, ls=":")
        m.plot([px], [.95], marker="v", color="#c0392b", ms=6); m.text(px, 1.08, p["pas_id"].split("_")[1], ha="center", fontsize=7, color="#c0392b")
    ts = truth_sites[r["tested_site_id"]]
    m.set_title(f"{g}  [{r['event_type']}]\n{ts['ctrl_frac']} -> {ts['kd_frac']} usage at {r['tested_site_id'].split('_')[1]}"
                + ("  (x3 expression)" if r["event_type"] == "expr_only" else ""), fontsize=9, loc="left", pad=12)
    ax[2, j].set_xlabel("position in gene, 5' -> 3' (bp)")
fig.suptitle("Synthetic IPA data: STAR-aligned coverage, mean of 4 replicates per condition (red markers = PAS reference)", fontsize=10, y=.995)
fig.tight_layout(); fig.savefig(out, dpi=140)
