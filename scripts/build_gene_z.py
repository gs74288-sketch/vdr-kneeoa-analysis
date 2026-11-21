# -*- coding: utf-8 -*-
"""Gene-level z-scores from knee OA GWAS summary statistics for scDRS.

Gene bodies are taken from a GTF; SNPs are assigned to any gene whose body
+/- 10 kb contains the SNP. Per gene the smallest p-value is kept and converted
to a signed z (sign from beta).

IMPORTANT: the GWAS file must use the same genome build as the GTF.
EBI GWAS Catalog harmonised files are GRCh38; using a GRCh37 GTF shifts genes by
up to ~1.5 Mb and silently produces wrong assignments (check a known locus such
as GDF5 before trusting the output).

Input : data/GCST007090.h.tsv.gz        (harmonised sumstats)
        data/Homo_sapiens.GRCh38.*.gtf.gz
Output: tables/GCST007090_gene_z.tsv    (columns: GENE, Z)
"""
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "tables")
os.makedirs(OUT, exist_ok=True)

GTF = os.path.join(DATA, "Homo_sapiens.GRCh38.112.gtf.gz")
GWAS = os.path.join(DATA, "GCST007090.h.tsv.gz")
OUTFILE = os.path.join(OUT, "GCST007090_gene_z.tsv")
PAD = 10000

print("parsing GTF ...", flush=True)
chrom_genes = {}
pat = re.compile(r'gene_id "([^".]+)(?:\.[0-9]+)?";.*?gene_name "([^"]+)";')
with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        chrom, start, end, attrs = parts[0], int(parts[3]), int(parts[4]), parts[8]
        m = pat.search(attrs)
        if m:
            gname = m.group(2)
        else:
            m2 = re.search(r'gene_id "([^"]+)"', attrs)
            if not m2:
                continue
            gname = m2.group(1).split(".")[0]
        chrom_genes.setdefault(chrom, []).append((start, end, gname))

for c in chrom_genes:
    chrom_genes[c].sort()
    arr = np.array([(s, e) for s, e, _ in chrom_genes[c]], dtype=np.int64)
    chrom_genes[c] = (arr[:, 0], arr[:, 1], [n for _, _, n in chrom_genes[c]])
print("genes:", sum(len(v[2]) for v in chrom_genes.values()), flush=True)

print("streaming GWAS ...", flush=True)
best = {}
nrow = 0
for chunk in pd.read_csv(GWAS, sep="\t", chunksize=500_000, low_memory=False):
    cols = {c.lower(): c for c in chunk.columns}
    chr_col = cols.get("hm_chrom") or cols.get("chromosome")
    pos_col = cols.get("hm_pos") or cols.get("base_pair_location")
    p_col = cols.get("p_value")
    b_col = cols.get("beta") or cols.get("hm_beta")
    if not (chr_col and pos_col and p_col):
        sys.exit("required columns not found: %s" % list(chunk.columns)[:10])
    d = chunk[[chr_col, pos_col, p_col] + ([b_col] if b_col else [])].copy()
    d.columns = ["chrom", "pos", "p"] + (["beta"] if b_col else [])
    d["p"] = pd.to_numeric(d["p"], errors="coerce")
    d["pos"] = pd.to_numeric(d["pos"], errors="coerce")
    d = d[d["p"].notna() & d["pos"].notna()]
    d["chrom"] = d["chrom"].astype(str).str.replace("chr", "", regex=False)
    nrow += len(d)
    for c, sub in d.groupby("chrom", sort=False):
        if c not in chrom_genes:
            continue
        starts, ends, names = chrom_genes[c]
        pos = sub["pos"].values.astype(np.int64)
        pv = sub["p"].values
        beta = sub["beta"].values if "beta" in sub else np.zeros(len(sub))
        j = np.searchsorted(starts, pos + PAD, side="right")
        for k in range(len(pos)):
            jj = j[k] - 1
            lo = pos[k] - PAD
            cnt = 0
            while jj >= 0 and cnt < 60:
                if ends[jj] >= lo:
                    g = names[jj]
                    cur = best.get(g)
                    if cur is None or pv[k] < cur[0]:
                        best[g] = (pv[k], beta[k])
                elif starts[jj] < lo - 1_000_000:
                    break
                jj -= 1
                cnt += 1
print("SNPs:", format(nrow, ","), "| genes with signal:", len(best), flush=True)

rows = []
for g, (p, b) in best.items():
    p = min(max(p, 1e-300), 1.0)
    z = stats.norm.isf(p / 2)
    if b is not None and b == b and b < 0:
        z = -z
    rows.append((g, z))
df = pd.DataFrame(rows, columns=["GENE", "Z"]).sort_values("Z", ascending=False)
df.to_csv(OUTFILE, sep="\t", index=False)
print("written:", OUTFILE, "| genes:", len(df), flush=True)

# sanity check on a known OA locus
for g in ("GDF5", "ALDH1A2"):
    r = df[df["GENE"] == g]
    if len(r):
        print("  %s: Z = %+.2f  (rank %d)" % (g, r["Z"].iloc[0],
                                              int((df["Z"].abs() > abs(r["Z"].iloc[0])).sum()) + 1))
