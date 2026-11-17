#!/usr/bin/env python3
"""Minimal reference caller used to sanity-check the synthetic data (NOT a replacement for your IPA tool).

For every polyA site: usage = site UMIs / all UMIs of that gene's (kept) sites, per replicate; Welch t-test on
logit(usage) between conditions; Benjamini-Hochberg across sites.  Call = |delta usage| >= --min-delta & padj < --alpha.
Intronic site up/down  -> IPA_gain / IPA_loss.   First (proximal) 3'UTR site of a 2-site gene up/down ->
3UTR_shortening / 3UTR_lengthening.   Sites whose downstream 20 nt hold >= --max-downA A's (internal priming
candidates) are dropped before usage is computed unless --no-ap-filter.
"""
import argparse, csv
import numpy as np
from scipy import stats

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="site_id<TAB>gene_id<TAB>sample... (count_site_umis.py output)")
    ap.add_argument("--pas", required=True); ap.add_argument("--meta", required=True)
    ap.add_argument("--ref", default="ctrl"); ap.add_argument("--target", default="KD")
    ap.add_argument("--min-delta", type=float, default=0.10); ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-umis", type=int, default=100, help="min total gene UMIs per sample to enter the test")
    ap.add_argument("--max-downA", type=int, default=12); ap.add_argument("--no-ap-filter", action="store_true")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    pas = {r["pas_id"]: r for r in csv.DictReader(open(a.pas), delimiter="\t")}
    cond = {r["sample"]: r["condition"] for r in csv.DictReader(open(a.meta))}
    rows = list(csv.reader(open(a.counts), delimiter="\t")); hdr, rows = rows[0], rows[1:]
    samples = hdr[2:]; ci = [i for i, s in enumerate(samples) if cond[s] == a.ref]; ti = [i for i, s in enumerate(samples) if cond[s] == a.target]
    by_gene = {}
    for r in rows:
        sid = r[0]
        if not a.no_ap_filter and int(pas[sid]["downstream_A20"]) >= a.max_downA: continue
        by_gene.setdefault(r[1], []).append((sid, np.array(r[2:], float)))
    res = []
    for g, sl in by_gene.items():
        tot = sum(c for _, c in sl)
        if len(sl) < 2 or (tot < a.min_umis).any(): continue
        for sid, c in sl:
            u = (c + 0.5) / (tot + 1.0)                      # pseudo-count keeps logit finite
            lg = np.log(u / (1 - u)); fr = c / tot
            p = stats.ttest_ind(lg[ti], lg[ci], equal_var=False).pvalue
            d = fr[ti].mean() - fr[ci].mean()
            res.append(dict(site_id=sid, gene_id=g, region=pas[sid]["region"], rank=int(pas[sid]["site_rank_5to3"]),
                            ctrl_frac=fr[ci].mean(), kd_frac=fr[ti].mean(), delta=d, pvalue=p, n_sites=len(sl)))
    p = np.array([r["pvalue"] for r in res]); o = np.argsort(p); m = len(p)
    q = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]; padj = np.empty(m); padj[o] = np.minimum(q, 1)
    for r, pa in zip(res, padj):
        r["padj"] = pa; call = "none"
        if pa < a.alpha and abs(r["delta"]) >= a.min_delta:
            up = r["delta"] > 0
            if r["region"] == "intronic": call = "IPA_gain" if up else "IPA_loss"
            elif r["n_sites"] >= 2 and pas[r["site_id"]]["site_rank_5to3"] != "" and r["site_id"].endswith("_utr1") \
                    and sum(1 for x in by_gene[r["gene_id"]] if x[0].split("_")[1].startswith("utr")) == 2:
                call = "3UTR_shortening" if up else "3UTR_lengthening"
        r["call"] = call
    with open(a.out, "w") as f:
        cols = ["site_id", "gene_id", "region", "ctrl_frac", "kd_frac", "delta", "pvalue", "padj", "call"]
        f.write("\t".join(cols) + "\n")
        for r in res: f.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"{len(res)} sites tested, {sum(r['call'] != 'none' for r in res)} called")

if __name__ == "__main__":
    main()
