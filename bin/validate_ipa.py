#!/usr/bin/env python3
"""Compare IPA / 3'UTR-APA calls with the planted truth.

  validate_ipa.py calls.tsv expected/truth_sites.tsv [--site-col site_id] [--call-col call]

calls.tsv needs a site id column (default `site_id`; matching `pas_id` is also accepted) and a call column whose
values are IPA_gain, IPA_loss, 3UTR_shortening, 3UTR_lengthening or none (rename with --call-col, or map your own
labels with --map "up=IPA_gain,down=IPA_loss").  A site missing from calls.tsv counts as `none`.

Reports: planted events found (with correct direction) per event type, and null / trap / expression-only sites
wrongly called.  Exit code 1 if sensitivity < --min-sens or false-positive rate > --max-fpr.
"""
import argparse, csv, sys
from collections import Counter, defaultdict

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("calls"); ap.add_argument("truth")
    ap.add_argument("--site-col", default=None); ap.add_argument("--call-col", default="call")
    ap.add_argument("--map", default="", help="label remap, e.g. 'up=IPA_gain,down=IPA_loss'")
    ap.add_argument("--min-sens", type=float, default=0.8); ap.add_argument("--max-fpr", type=float, default=0.10)
    a = ap.parse_args()
    rm = dict(x.split("=", 1) for x in a.map.split(",") if x)
    sep = "\t" if a.calls.endswith((".tsv", ".txt")) else ","
    cr = csv.DictReader(open(a.calls), delimiter=sep)
    sc = a.site_col or next((c for c in ("site_id", "pas_id") if c in cr.fieldnames), None)
    if sc is None or a.call_col not in cr.fieldnames: sys.exit(f"need site + '{a.call_col}' columns; found {cr.fieldnames}")
    calls = {r[sc]: rm.get(r[a.call_col], r[a.call_col]) or "none" for r in cr}
    truth = list(csv.DictReader(open(a.truth), delimiter="\t"))
    ev = defaultdict(Counter); null = defaultdict(Counter); tp = pos = fp = neg = 0
    for t in truth:
        got = calls.get(t["site_id"], "none"); exp = t["expected_call"]
        if exp != "none":
            pos += 1; ok = got == exp; tp += ok
            ev[t["event_type"]]["found" if ok else ("wrong_direction" if got != "none" else "missed")] += 1
        elif t["site_role"] in ("tested_site", "artifact_site"):           # one row per null gene / trap
            neg += 1; bad = got != "none"; fp += bad
            null[t["event_type"]]["false_call" if bad else "clean"] += 1
    print(f"{'event type':<12}{'found':>7}{'wrong dir':>11}{'missed':>8}")
    for k, c in sorted(ev.items()): print(f"{k:<12}{c['found']:>7}{c['wrong_direction']:>11}{c['missed']:>8}")
    print(f"\n{'null type':<12}{'clean':>7}{'false call':>12}")
    for k, c in sorted(null.items()): print(f"{k:<12}{c['clean']:>7}{c['false_call']:>12}")
    sens, fpr = tp / max(pos, 1), fp / max(neg, 1)
    print(f"\nplanted events found: {tp}/{pos} ({sens:.0%});  null/trap/expr-only genes wrongly called: {fp}/{neg} ({fpr:.0%})")
    sys.exit(0 if sens >= a.min_sens and fpr <= a.max_fpr else 1)

if __name__ == "__main__":
    main()
