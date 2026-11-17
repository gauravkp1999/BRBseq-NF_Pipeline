#!/usr/bin/env python3
"""Compare the pipeline's merged count matrix with the synthetic ground truth.

usage: validate_counts.py results/counts/brbseq.merged_counts.tsv assets/test_data/expected/expected_counts.tsv
Exit code 0 = PASS, 1 = FAIL.
"""
import sys
import numpy as np
import pandas as pd

obs = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
exp = pd.read_csv(sys.argv[2], sep="\t", index_col=0)
ok = True

missing_cols = sorted(set(exp.columns) - set(obs.columns))
extra_cols = sorted(set(obs.columns) - set(exp.columns))
if missing_cols or extra_cols:
    print(f"FAIL sample columns differ. missing={missing_cols} unexpected={extra_cols}")
    ok = False
common_cols = [c for c in exp.columns if c in obs.columns]
obs = obs.reindex(exp.index).fillna(0)[common_cols]
exp = exp[common_cols]

print(f"genes: expected {len(exp)}, matrix rows matched {int((obs.sum(axis=1) > 0).sum())}")
print(f"{'sample':8}{'truth_total':>12}{'observed':>10}{'ratio':>8}{'pearson':>9}{'exact_genes':>13}  best_match")
for c in common_cols:
    r = np.corrcoef(obs[c], exp[c])[0, 1] if obs[c].std() and exp[c].std() else float("nan")
    exact = float((obs[c] == exp[c]).mean())
    ratio = obs[c].sum() / max(exp[c].sum(), 1)
    # detect sample swaps: which truth column correlates best with this observed column?
    best = max(common_cols, key=lambda e: np.corrcoef(obs[c], exp[e])[0, 1] if obs[c].std() else -1)
    flag = "" if best == c else "  <-- SWAPPED?"
    print(f"{c:8}{int(exp[c].sum()):>12}{int(obs[c].sum()):>10}{ratio:>8.3f}{r:>9.4f}{exact:>13.3f}  {best}{flag}")
    if best != c or r < 0.99 or not (0.97 <= ratio <= 1.03):
        ok = False

print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
