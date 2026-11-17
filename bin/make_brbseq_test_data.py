#!/usr/bin/env python3
"""Generate a SYNTHETIC BRB-seq dataset for testing intronic polyadenylation (IPA) and 3'UTR-APA analysis.

Input format is identical to the gene-level brbseq pipeline test set (R1 = 12 bp barcode + 16 bp UMI,
R2 = 80 bp cDNA, one barcode map per run, samplesheet.csv), so it can be fed to the same STARsolo front end.

Reference
---------
* 3 chromosomes, ~92 multi-exon genes (3-5 exons, 1.1-1.7 kb introns, GT..AG splice motifs), both strands.
* Every gene has a terminal exon carrying the normal polyA site(s); most also carry >=1 intronic polyA site.
* PAS reference (TSV + BED): every candidate site the tool should test, with hexamer / downstream-A columns.

Reads
-----
* 3'-biased cDNA, sense strand, ending 0-250 nt upstream of the cleavage site (~12% of reads run into the
  non-templated poly(A) tail).  Intronic-site reads are unspliced and end inside the intron.
* ~3% of each gene's UMIs are 3'-bias-free background reads from the spliced transcript.
* Internal-priming trap genes: an A-rich genomic stretch (>=17/20 A) inside an intron; reads end just
  upstream of it and carry no non-templated tail.  They are NOT real polyA sites.

Planted design (fractions = share of the gene's UMIs using the site)
---------------------------------------------------------------------
  ipa_up        intronic  ~5%  -> ~40%   ctrl -> KD
  ipa_down      intronic  ~40% -> ~5%
  ipa_null      intronic  8-30%, identical in both conditions (negative control)
  utr_short     proximal 3'UTR site ~17% -> ~70%  (3'UTR shortening, no intronic site)
  utr_long      proximal 3'UTR site ~70% -> ~17%  (3'UTR lengthening)
  utr_null      proximal 3'UTR site 30-50%, identical
  expr_only     gene x3 in KD, site fractions identical (must NOT be called IPA)
  trap          internal-priming A-stretch; half of the genes use it more in KD (false IPA if unfiltered)
Replicate noise: library size (lognormal), NB gene counts, lognormal noise on every site weight.

Usage: python make_ipa_test_data.py --outdir assets/test_data_ipa [--reps 4] [--seed 21] [--mean-umis 550]
"""
import argparse, gzip, json, os, random
import numpy as np

COMP = str.maketrans("ACGT", "TGCA")
LOW_A = [0.15, 0.285, 0.285, 0.28]


def rc(s):
    return s.translate(COMP)[::-1]


def rseq(rng, n, w=None):
    return "".join(rng.choices("ACGT", weights=w, k=n))


def ham(a, b):
    return sum(x != y for x, y in zip(a, b))


def barcodes(rng, n, length, min_dist):
    out = []
    while len(out) < n:
        c = rseq(rng, length)
        if all(ham(c, o) >= min_dist for o in out):
            out.append(c)
    return out


def write_gz(path, text):
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as f:
        f.write(text.encode())


def neighbours1(u):
    for i in range(len(u)):
        for b in "ACGT":
            if b != u[i]:
                yield u[:i] + b + u[i + 1:]


def unique_umis(rng, n, length=16):
    out, seen = [], set()
    while len(out) < n:
        u = rseq(rng, length)
        if u in seen or any(v in seen for v in neighbours1(u)):
            continue
        seen.add(u)
        out.append(u)
    return out


# --------------------------------------------------------------------------------------- categories
CATS = [("ipa_up", 12), ("ipa_down", 12), ("ipa_null", 12), ("utr_short", 12), ("utr_long", 12),
        ("utr_null", 12), ("expr_only", 12), ("trap", 8)]
UTR_CATS = {"utr_short", "utr_long", "utr_null"}
CALL = {"ipa_up": "IPA_gain", "ipa_down": "IPA_loss", "utr_short": "3UTR_shortening", "utr_long": "3UTR_lengthening"}


def fractions(rng, cat):
    """Per-condition target fractions for the *changing / tested* site. Returns (ctrl, kd, trap_kind)."""
    u = rng.uniform
    if cat == "ipa_up":    return u(.03, .08), u(.32, .48), None
    if cat == "ipa_down":  return u(.32, .48), u(.03, .08), None
    if cat == "ipa_null":  f = u(.08, .30); return f, f, None
    if cat == "expr_only": f = u(.10, .30); return f, f, None
    if cat == "utr_short": return u(.10, .25), u(.60, .80), None
    if cat == "utr_long":  return u(.60, .80), u(.10, .25), None
    if cat == "utr_null":  f = u(.30, .50); return f, f, None
    raise ValueError(cat)


def build_gene(rng, gid, cat, strand):
    n_ex = rng.choice([3, 4, 5])
    utr = cat in UTR_CATS
    ex_len = [rng.randint(180, 320)] + [rng.randint(110, 200) for _ in range(n_ex - 2)] + \
             [1300 if utr else rng.randint(800, 1000)]
    in_len = [rng.randint(1100, 1700) for _ in range(n_ex - 1)]
    parts, ex, intr, pos = [], [], [], 0
    for i in range(n_ex):
        parts.append(rseq(rng, ex_len[i]))
        ex.append((pos, pos + ex_len[i])); pos += ex_len[i]
        if i < n_ex - 1:
            L = in_len[i]
            parts.append("GT" + rseq(rng, L - 24) + rseq(rng, 20, [.05, .45, .05, .45]) + "AG")
            intr.append((pos, pos + L)); pos += L
    flank = 60
    parts.append(rseq(rng, flank, LOW_A))
    seq = list("".join(parts))
    ts, te = ex[-1]

    def intron_pick(exclude=None):
        w = [0.2] * (n_ex - 1); w[-1] = 0.55
        while True:
            k = rng.choices(range(n_ex - 1), weights=w)[0]
            if exclude is None or k != exclude:
                return k

    sites = []   # dict(id, kind, off, k)  (off = gene-oriented 0-based offset of last transcribed base)
    if utr:
        sites.append(dict(kind="utr", off=ts + rng.randint(450, 600)))     # proximal
        sites.append(dict(kind="utr", off=te - 1))                         # distal
    else:
        sites.append(dict(kind="utr", off=te - 1))
    n_int, trap = 0, None
    if cat in ("ipa_up", "ipa_down", "ipa_null", "expr_only"):
        k = intron_pick(); s0, e0 = intr[k]
        sites.append(dict(kind="intronic", off=s0 + rng.randint(450, (e0 - s0) - 250), k=k, changing=True))
        if cat != "expr_only" and rng.random() < 0.35:                    # extra constant intronic site
            for _ in range(50):
                k2 = intron_pick(); s2, e2 = intr[k2]
                o2 = s2 + rng.randint(450, (e2 - s2) - 250)
                if abs(o2 - sites[-1]["off"]) > 300 or k2 != k:
                    sites.append(dict(kind="intronic", off=o2, k=k2, changing=False)); break
    if cat == "trap":
        k = intron_pick(); s0, e0 = intr[k]
        off = s0 + rng.randint(450, (e0 - s0) - 300)
        sites.append(dict(kind="trap", off=off, k=k, changing=True))
    # signals around real sites
    for s in sites:
        P = s["off"]
        if s["kind"] == "trap":
            st = list("AAAAAAAAAAAAAAAAAAAA")
            for j in rng.sample(range(20), 2): st[j] = rng.choice("CGT")
            seq[P + 1:P + 21] = st
            s["hexamer"] = "NA"
            continue
        hx = "AATAAA" if rng.random() < 0.85 else "ATTAAA"
        seq[P - 24:P - 18] = list(hx); s["hexamer"] = hx
        while seq[P + 1:P + 21].count("A") > 6:
            seq[P + 1:P + 21] = list(rseq(rng, 20, LOW_A))
    seq = "".join(seq)
    for s in sites: s["downA20"] = seq[s["off"] + 1:s["off"] + 21].count("A")
    # ordering / ids
    utr_sites = sorted([s for s in sites if s["kind"] == "utr"], key=lambda s: s["off"])
    for i, s in enumerate(utr_sites, 1): s["id"] = f"{gid}_utr{i}"
    int_sites = sorted([s for s in sites if s["kind"] != "utr"], key=lambda s: s["off"])
    for i, s in enumerate(int_sites, 1): s["id"] = f"{gid}_int{i}"
    # fractions
    g = dict(id=gid, cat=cat, strand=strand, seq=seq, exons=ex, introns=intr, sites=sites, flank=flank,
             stop_off=ts + 100)
    g["frac"] = {"ctrl": {}, "KD": {}}
    if cat == "trap":
        diff = None  # decided by caller (alternating) via g["trap_diff"]
    else:
        fc, fk, _ = fractions(rng, cat)
        chg = [s for s in sites if (s["kind"] == "utr" and utr and s is utr_sites[0]) or
               (s["kind"] == "intronic" and s.get("changing"))][0]
        const_int = [s for s in sites if s["kind"] == "intronic" and not s.get("changing")]
        const_f = {s["id"]: rng.uniform(.04, .08) for s in const_int}
        for cond, f in (("ctrl", fc), ("KD", fk)):
            d = {chg["id"]: f, **const_f}
            d[utr_sites[-1]["id"]] = 1.0 - sum(d.values())
            g["frac"][cond] = d
        g["changing_id"] = chg["id"]
    return g


def tx_seq(g, site=None):
    """Spliced (gene-oriented) transcript ending at the site; site=None -> full-length mature transcript."""
    seq, ex = g["seq"], g["exons"]
    if site is None:
        return "".join(seq[s:e] for s, e in ex)
    if site["kind"] == "utr":
        return "".join(seq[s:e] for s, e in ex[:-1]) + seq[ex[-1][0]:site["off"] + 1]
    k = site["k"]
    return "".join(seq[s:e] for s, e in ex[:k + 1]) + seq[g["introns"][k][0]:site["off"] + 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="assets/test_data_ipa")
    ap.add_argument("--rel-prefix", default=None, help="path prefix written into samplesheet.csv (default: --outdir)")
    ap.add_argument("--reps", type=int, default=4, help="replicates per condition")
    ap.add_argument("--seed", type=int, default=21)
    ap.add_argument("--mean-umis", type=float, default=550.0, help="median unique UMIs per gene per sample")
    ap.add_argument("--polya-frac", type=float, default=0.12, help="fraction of site reads that run into the poly(A) tail")
    ap.add_argument("--site-sd", type=float, default=0.30, help="replicate noise (log-scale) on each site weight")
    a = ap.parse_args()
    rng = random.Random(a.seed); nprng = np.random.default_rng(a.seed)
    out, prefix = a.outdir, a.rel_prefix or a.outdir
    for d in ("reference", "expected"): os.makedirs(f"{out}/{d}", exist_ok=True)

    # ------------------------------------------------------------------ genes
    cats = [c for c, n in CATS for _ in range(n)]
    rng.shuffle(cats)
    genes, ti = [], 0
    for i, cat in enumerate(cats):
        g = build_gene(rng, f"GENE{i + 1:04d}", cat, rng.choice("++-"))
        if cat == "trap":
            g["trap_diff"] = (ti % 2 == 0); ti += 1
            tr = [s for s in g["sites"] if s["kind"] == "trap"][0]
            utr1 = [s for s in g["sites"] if s["kind"] == "utr"][0]
            if g["trap_diff"]: fc, fk = rng.uniform(.02, .05), rng.uniform(.20, .30)
            else:              fc = fk = rng.uniform(.10, .15)
            for cond, f in (("ctrl", fc), ("KD", fk)):
                g["frac"][cond] = {tr["id"]: f, utr1["id"]: 1 - f}
            g["changing_id"] = tr["id"]
        genes.append(g)
    nchr = 3
    for i, g in enumerate(genes): g["chrom"] = f"chr{i * nchr // len(genes) + 1}"

    # ------------------------------------------------------------------ genome + GTF + PAS reference
    chroms, gtf, pas = {}, [], []
    for chrom in sorted({g["chrom"] for g in genes}):
        parts, pos = [rseq(rng, 400)], 400
        for g in [g for g in genes if g["chrom"] == chrom]:
            block = g["seq"]; Lb = len(block)
            g["G0"] = pos; g["Lb"] = Lb
            parts.append(block if g["strand"] == "+" else rc(block))
            sp = rseq(rng, rng.randint(350, 550)); parts.append(sp); pos += Lb + len(sp)
        chroms[chrom] = "".join(parts)
    for g in genes:
        G0, Lb = g["G0"], g["Lb"]
        if g["strand"] == "+":
            iv = lambda s, e: (G0 + s + 1, G0 + e)
            gpos = lambda o: G0 + o + 1
        else:
            iv = lambda s, e: (G0 + Lb - e + 1, G0 + Lb - s)
            gpos = lambda o: G0 + Lb - o
        g["iv"], g["gpos"] = iv, gpos
        ex = g["exons"]
        a0, a1 = iv(ex[0][0], ex[-1][1])
        attr = f'gene_id "{g["id"]}"; transcript_id "{g["id"]}.1"; gene_name "{g["id"]}"; gene_biotype "protein_coding";'
        gtf.append((g["chrom"], a0, "gene", a1, g["strand"], attr))
        gtf.append((g["chrom"], a0, "transcript", a1, g["strand"], attr))
        for n, (s, e) in enumerate(ex, 1):
            x0, x1 = iv(s, e)
            gtf.append((g["chrom"], x0, "exon", x1, g["strand"], attr + f' exon_number "{n}";'))
        u0, u1 = iv(g["stop_off"] + 1, ex[-1][1])
        gtf.append((g["chrom"], u0, "three_prime_utr", u1, g["strand"], attr))
        # sanity: GTF exon coordinates index the right bases
        for (s, e) in ex:
            x0, x1 = iv(s, e); sub = chroms[g["chrom"]][x0 - 1:x1]
            assert (sub if g["strand"] == "+" else rc(sub)) == g["seq"][s:e]
        g["sites"].sort(key=lambda s: s["off"])
        for rank, s in enumerate(g["sites"], 1):
            s["pos"] = gpos(s["off"])
            if s["kind"] == "utr":
                reg, loc = "3UTR", f"exon{len(ex)}"
            else:
                reg, loc = "intronic", f"intron{s['k'] + 1}"
            s["region"], s["loc"] = reg, loc
            pas.append(dict(pas_id=s["id"], gene_id=g["id"], gene_name=g["id"], chrom=g["chrom"], strand=g["strand"],
                            pos=s["pos"], region=reg, feature=loc, is_last_exon=int(s["kind"] == "utr"),
                            site_rank_5to3=rank, n_sites_in_gene=len(g["sites"]), hexamer=s["hexamer"],
                            downstream_A20=s["downA20"]))
    gtf.sort(key=lambda r: (r[0], r[1], {"gene": 0, "transcript": 1}.get(r[2], 2)))
    with open(f"{out}/reference/genome.fasta", "w") as f:
        for c, s in chroms.items():
            f.write(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    with open(f"{out}/reference/genes.gtf", "w") as f:
        for c, st, ft, en, sd, at in gtf: f.write(f"{c}\tsynthetic\t{ft}\t{st}\t{en}\t.\t{sd}\t.\t{at}\n")
    cols = list(pas[0].keys())
    with open(f"{out}/reference/pas_reference.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in pas: f.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(f"{out}/reference/pas_reference.bed", "w") as f:
        for r in sorted(pas, key=lambda r: (r["chrom"], r["pos"])):
            f.write(f"{r['chrom']}\t{r['pos'] - 1}\t{r['pos']}\t{r['pas_id']}\t{r['downstream_A20']}\t{r['strand']}\n")

    # ------------------------------------------------------------------ design
    conds = {"ctrl": [f"ctrl_{i}" for i in range(1, a.reps + 1)], "KD": [f"KD_{i}" for i in range(1, a.reps + 1)]}
    samples = conds["ctrl"] + conds["KD"]
    cond_of = {s: c for c, ss in conds.items() for s in ss}
    plates = {"plate1": [], "plate2": []}
    for c in ("ctrl", "KD"):
        for i, s in enumerate(conds[c]): plates[f"plate{i % 2 + 1}"].append(s)
    bcs = barcodes(rng, max(len(v) for v in plates.values()), 12, 4)
    size = {s: float(np.exp(nprng.normal(0, 0.18))) for s in samples}
    for g in genes:
        g["mu"] = float(np.clip(np.exp(nprng.normal(np.log(a.mean_umis), 0.45)), 120, 3000))

    # ------------------------------------------------------------------ simulate molecules
    BG = 0.03
    site_counts = {}                    # (site_id, sample) -> unique UMIs
    gene_total = {}
    mols = {p: [] for p in plates}      # (sample, r2_seq, umi, n_copies)
    for s in samples:
        plate = next(p for p, ss in plates.items() if s in ss)
        for g in genes:
            c = cond_of[s]
            mu = g["mu"] * size[s] * (3.0 if (c == "KD" and g["cat"] == "expr_only") else 1.0)
            phi = 0.05
            total = int(nprng.poisson(nprng.gamma(1 / phi, mu * phi)))
            ids = list(g["frac"][c].keys())
            w = np.array([g["frac"][c][i] for i in ids]) * np.exp(nprng.normal(0, a.site_sd, len(ids)))
            p = np.append(w / w.sum() * (1 - BG), BG)
            cnt = nprng.multinomial(total, p)
            gene_total[(g["id"], s)] = total
            umis = unique_umis(rng, total)
            idx = 0
            txs = {sid: tx_seq(g, next(x for x in g["sites"] if x["id"] == sid)) for sid in ids}
            full = tx_seq(g)
            for sid, n in zip(ids + [g["id"] + "_bg"], cnt):
                site_counts[(sid, s)] = int(n)
                for _ in range(int(n)):
                    u = umis[idx]; idx += 1
                    if sid.endswith("_bg"):
                        L = len(full); st = rng.randint(0, L - 80); r2 = full[st:st + 80]
                    else:
                        tx = txs[sid]; L = len(tx)
                        if rng.random() < a.polya_frac:                 # read runs into the poly(A) tail
                            ext = rng.randint(4, 35); r2 = (tx + "A" * 40)[L + ext - 80:L + ext]
                        else:
                            d = min(int(nprng.exponential(45)), 250); r2 = tx[L - d - 80:L - d]
                    if g["strand"] == "-": pass                           # tx is already in sense orientation
                    copies = min(4, int(nprng.geometric(0.7)))
                    mols[plate].append((s, r2, u, copies))

    # ------------------------------------------------------------------ write reads
    ERR = 0.001
    stats = {"total_reads": 0, "valid_barcode_reads": 0, "invalid_barcode_reads": 0, "mismatch_barcode_reads": 0}
    for plate, ss in plates.items():
        bc_of = {s: bcs[i] for i, s in enumerate(ss)}
        recs = []
        for s, r2, u, copies in mols[plate]:
            for _ in range(copies):
                b = bc_of[s]
                if rng.random() < 0.015:
                    i = rng.randrange(12); b = b[:i] + rng.choice([x for x in "ACGT" if x != b[i]]) + b[i + 1:]
                    stats["mismatch_barcode_reads"] += 1
                q = ["I"] * 80; seq = list(r2)
                for j in range(80):
                    if rng.random() < ERR:
                        seq[j] = rng.choice([x for x in "ACGT" if x != seq[j]]); q[j] = "5"
                recs.append((b + u, "".join(seq), "".join(q)))
                stats["valid_barcode_reads"] += 1
        for _ in range(int(len(recs) * 0.02)):
            while True:
                b = rseq(rng, 12)
                if all(ham(b, x) >= 2 for x in bcs): break
            g = rng.choice(genes); tx = tx_seq(g); st = rng.randint(0, len(tx) - 80)
            recs.append((b + rseq(rng, 16), tx[st:st + 80], "I" * 80)); stats["invalid_barcode_reads"] += 1
        rng.shuffle(recs); stats["total_reads"] += len(recs)
        write_gz(f"{out}/{plate}_R1.fastq.gz", "".join(f"@{plate}_{i}\n{x}\n+\n{'I' * len(x)}\n" for i, (x, _, _) in enumerate(recs)))
        write_gz(f"{out}/{plate}_R2.fastq.gz", "".join(f"@{plate}_{i}\n{y}\n+\n{q}\n" for i, (_, y, q) in enumerate(recs)))
        with open(f"{out}/barcode_map_{plate}.csv", "w") as f:
            f.write("sample,barcode\n" + "".join(f"{s},{bc_of[s]}\n" for s in ss))

    # ------------------------------------------------------------------ sheets
    with open(f"{out}/samplesheet.csv", "w") as f:
        f.write("run,fastq_1,fastq_2,barcode_map\n")
        for p in plates:
            f.write(f"{p},{prefix}/{p}_R1.fastq.gz,{prefix}/{p}_R2.fastq.gz,{prefix}/barcode_map_{p}.csv\n")
    with open(f"{out}/sample_meta.csv", "w") as f:
        f.write("sample,condition,plate\n")
        for p, ss in plates.items():
            for s in ss: f.write(f"{s},{cond_of[s]},{p}\n")
    with open(f"{out}/contrasts.csv", "w") as f:
        f.write("id,variable,reference,target\nKD_vs_ctrl,condition,ctrl,KD\n")
    ids = [g["id"] for g in genes]; eo = [g["id"] for g in genes if g["cat"] == "expr_only"]
    with open(f"{out}/genesets.gmt", "w") as f:
        f.write("UP_IN_KD\tgenes planted x3 in KD\t" + "\t".join(eo + rng.sample([i for i in ids if i not in eo], 5)) + "\n")
        for i in range(1, 7):
            f.write(f"RANDOM_SET_{i}\trandom genes\t" + "\t".join(rng.sample(ids, rng.randint(15, 40))) + "\n")

    # ------------------------------------------------------------------ truth
    order = {}
    sites_all = [(g, s) for g in genes for s in g["sites"]]
    with open(f"{out}/expected/expected_site_counts.tsv", "w") as f:
        f.write("site_id\tgene_id\t" + "\t".join(samples) + "\n")
        for g, s in sites_all:
            f.write(f"{s['id']}\t{g['id']}\t" + "\t".join(str(site_counts[(s['id'], x)]) for x in samples) + "\n")
        for g in genes:
            f.write(f"{g['id']}_bg\t{g['id']}\t" + "\t".join(str(site_counts[(g['id'] + '_bg', x)]) for x in samples) + "\n")
    with open(f"{out}/expected/expected_gene_counts_exonic.tsv", "w") as f:       # what STARsolo 'Gene' can see
        f.write("gene_id\t" + "\t".join(samples) + "\n")
        for g in genes:
            keep = [s["id"] for s in g["sites"] if s["kind"] == "utr"] + [g["id"] + "_bg"]
            f.write(g["id"] + "\t" + "\t".join(str(sum(site_counts[(k, x)] for k in keep)) for x in samples) + "\n")
    with open(f"{out}/expected/expected_gene_counts_total.tsv", "w") as f:
        f.write("gene_id\t" + "\t".join(samples) + "\n")
        for g in genes: f.write(g["id"] + "\t" + "\t".join(str(gene_total[(g["id"], x)]) for x in samples) + "\n")

    hdr = ["gene_id", "site_id", "chrom", "strand", "position", "region", "event_type", "site_role",
           "ctrl_frac", "kd_frac", "delta_frac", "expected_call", "is_null", "gene_expr_fc_KD_vs_ctrl", "note"]
    notes = {"ipa_up": "intronic site gains usage in KD", "ipa_down": "intronic site loses usage in KD",
             "ipa_null": "intronic site, unchanged: negative control",
             "utr_short": "proximal 3'UTR site gains usage in KD", "utr_long": "proximal 3'UTR site loses usage in KD",
             "utr_null": "two 3'UTR sites, unchanged: negative control",
             "expr_only": "whole gene x3 in KD, site ratios unchanged: must not be called IPA"}
    with open(f"{out}/expected/truth_sites.tsv", "w") as f, open(f"{out}/expected/truth_genes.tsv", "w") as fg:
        f.write("\t".join(hdr) + "\n")
        fg.write("gene_id\tevent_type\texpected_call\tis_null\tgene_expr_fc_KD_vs_ctrl\ttested_site_id\n")
        for g in genes:
            cat = g["cat"]; fcx = 3 if cat == "expr_only" else 1
            if cat == "trap":
                call = "none"; null = True
                note = ("internal-priming artifact: usage rises in KD, so a missing A-filter gives a FALSE IPA gain"
                        if g["trap_diff"] else "internal-priming artifact, constant usage; must be filtered")
            else:
                call = CALL.get(cat, "none"); null = cat not in CALL; note = notes[cat]
            fg.write(f"{g['id']}\t{cat}\t{call}\t{str(null).upper()}\t{fcx}\t{g['changing_id']}\n")
            for s in g["sites"]:
                c0, k0 = g["frac"]["ctrl"][s["id"]], g["frac"]["KD"][s["id"]]
                if s["id"] == g["changing_id"]:
                    role = "artifact_site" if cat == "trap" else "tested_site"; sc = call
                    if cat == "trap" and s["kind"] == "trap": sc = "none"
                else:
                    role = "constant_site" if abs(c0 - k0) < 1e-9 else "compensating_site"; sc = "none"
                    if role == "compensating_site" and cat in CALL: role = "compensating_site"
                f.write("\t".join(map(str, [g["id"], s["id"], g["chrom"], g["strand"], s["pos"], s["region"], cat, role,
                                            f"{c0:.4f}", f"{k0:.4f}", f"{k0 - c0:+.4f}", sc,
                                            str(sc == "none").upper(), fcx, note])) + "\n")
    cnt = {c: sum(1 for g in genes if g["cat"] == c) for c, _ in CATS}
    stats.update(seed=a.seed, n_genes=len(genes), n_sites=len(sites_all), reps=a.reps, samples=samples,
                 genes_per_category=cnt, genome_bp=sum(len(s) for s in chroms.values()),
                 expected_valid_barcode_fraction=round(stats["valid_barcode_reads"] / stats["total_reads"], 4))
    json.dump(stats, open(f"{out}/expected/truth.json", "w"), indent=2)
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
