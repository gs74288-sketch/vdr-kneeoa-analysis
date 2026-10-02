# -*- coding: utf-8 -*-
"""Directional replication of VDR detection enrichment in regulatory/homeostatic
chondrocyte states, independent dataset GSE324993 (OA + healthy cartilage scRNA-seq).
Handles a gzip stream whose final trailer is missing at source (zlib streaming)."""
import gzip, csv, json, math
import zlib
import numpy as np
import pandas as pd
from scipy import stats

BASE = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002"
OUT = BASE + r"\outputs"
MAT = BASE + r"\data\GSE324993_counts_matrix.csv.gz"
ANN = BASE + r"\data\GSE324993_cell_annotation.txt.gz"
META = BASE + r"\data\GSE324993_sample_metadata.txt.gz"
GENES = ["VDR", "COL2A1", "ACAN", "COL10A1", "MMP13", "COL1A1", "SOX9", "PRG4", "IBSP", "MKI67"]

# annotation
lab2type = {}
with gzip.open(ANN, "rt", encoding="utf-8", errors="replace") as f:
    started = False
    for line in f:
        line = line.rstrip("\n")
        if not started:
            if line.startswith("Cell_Label"):
                started = True
            continue
        if not line.strip():
            continue
        p = line.split("\t")
        if len(p) >= 2:
            lab2type[p[0]] = p[1]

samp2group = {}
with gzip.open(META, "rt", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    next(rd)
    for row in rd:
        if len(row) >= 2:
            samp2group[row[0]] = row[1]

print("annotated cells:", len(lab2type), "samples:", samp2group, flush=True)

# ---- stream-decompress with zlib (source gz lacks trailer) ----
dec = zlib.decompressobj(31)
cells = None
totals = None
gene_rows = {}
buf = ""
nrow = 0
with open(MAT, "rb") as fh:
    while True:
        chunk = fh.read(1 << 20)
        if not chunk:
            break
        try:
            buf += dec.decompress(chunk).decode("utf-8", errors="replace")
        except zlib.error as e:
            print("zlib stopped:", e, flush=True)
            break
        *lines, buf = buf.split("\n")
        for line in lines:
            line = line.rstrip("\r")
            if not line:
                continue
            if cells is None:
                cells = line.split(",")
                totals = np.zeros(len(cells), dtype=np.int64)
                print("cells:", len(cells), flush=True)
                continue
            i = line.find(",")
            if i < 0:
                continue
            g = line[:i]
            vals = np.fromiter(map(int, line[i+1:].split(",")), dtype=np.int64, count=len(cells))
            if len(vals) != len(cells):
                continue
            totals += vals
            if g in GENES:
                gene_rows[g] = vals.astype(np.float32)
            nrow += 1
            if nrow % 5000 == 0:
                print("genes parsed:", nrow, flush=True)

print("total genes parsed:", nrow, "captured:", list(gene_rows.keys()), flush=True)

vdr = gene_rows.get("VDR")
df = pd.DataFrame({"cell": cells, "vdr": vdr, "total": totals})
df["sample"] = df["cell"].map(lambda c: c.split("_")[0])
df["group"] = df["sample"].map(lambda s: samp2group.get(s, "NA"))
df["cell_type"] = df["cell"].map(lambda c: lab2type.get(c, "NA"))
df["detected"] = (df["vdr"] > 0).astype(int)

def base_state(ct):
    if ct == "NA":
        return "NA"
    s = ct.split("(")[0].strip()
    for st in ["RegC", "HomC", "HTC", "preHTC", "ProC", "EC", "FC", "MTC"]:
        if s.startswith(st):
            return st
    return s

df["state"] = df["cell_type"].map(base_state)

summary = df.groupby("state").agg(n=("cell", "size"), detected=("detected", "sum"),
                                  mean_total=("total", "mean")).reset_index()
summary["frac"] = summary["detected"] / summary["n"]
summary = summary.sort_values("frac", ascending=False)
summary.to_csv(OUT + r"\GSE324993_VDR_by_state.csv", index=False)
print("\nVDR detection fraction by state:\n", summary.to_string(index=False), flush=True)

states = ["RegC", "HomC", "HTC", "preHTC", "ProC", "EC", "FC", "MTC"]
res = {}
for grp, label in [("osteoarthritis", "OA"), ("healthy controls", "HC"), (None, "ALL")]:
    sub = df if grp is None else df[df["group"] == grp]
    sub = sub[sub["state"].isin(states)]
    reg = sub[sub["state"].isin(["RegC", "HomC"])]
    oth = sub[~sub["state"].isin(["RegC", "HomC"])]
    if len(reg) and len(oth):
        table = [[int(reg["detected"].sum()), int(len(reg)-reg["detected"].sum())],
                 [int(oth["detected"].sum()), int(len(oth)-oth["detected"].sum())]]
        orr, p = stats.fisher_exact(table)
        res[label] = dict(reg_n=int(len(reg)), reg_frac=float(reg["detected"].mean()),
                          other_n=int(len(oth)), other_frac=float(oth["detected"].mean()),
                          OR=float(orr), P=float(p))
        print(f"[{label}] RegC/HomC {len(reg)} ({reg['detected'].mean():.4f}) vs Other "
              f"{len(oth)} ({oth['detected'].mean():.4f}); OR={orr:.3f} P={p:.4g}", flush=True)

per_donor = []
for s in sorted(df["sample"].unique()):
    sub = df[df["sample"] == s]
    reg = sub[sub["state"].isin(["RegC", "HomC"])]
    oth = sub[~sub["state"].isin(["RegC", "HomC"])]
    if len(reg) and len(oth):
        per_donor.append(dict(sample=s, group=sub["group"].iloc[0], reg_n=int(len(reg)),
                              reg_frac=float(reg["detected"].mean()), other_n=int(len(oth)),
                              other_frac=float(oth["detected"].mean()),
                              delta=float(reg["detected"].mean()-oth["detected"].mean())))
pd.DataFrame(per_donor).to_csv(OUT + r"\GSE324993_VDR_per_donor.csv", index=False)
print("\nPer-donor deltas:", flush=True)
for r in per_donor:
    print("  ", r, flush=True)

# depth-adjusted logistic (IRLS)
def logistic(X, y, iters=60):
    beta = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = X @ beta; mu = 1/(1+np.exp(-eta))
        w = np.clip(mu*(1-mu), 1e-9, None)
        z = eta + (y-mu)/w
        XtW = X.T * w
        try:
            bn = np.linalg.solve(XtW @ X, XtW @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(bn-beta)) < 1e-9:
            beta = bn; break
        beta = bn
    eta = X @ beta; mu = 1/(1+np.exp(-eta)); w = np.clip(mu*(1-mu), 1e-9, None)
    cov = np.linalg.inv((X*w[:, None]).T @ X)
    return beta, np.sqrt(np.diag(cov))

sub = df[df["state"].isin(states)].copy()
sub["is_reg"] = sub["state"].isin(["RegC", "HomC"]).astype(float)
sub["log_total"] = np.log1p(sub["total"])
X = np.column_stack([np.ones(len(sub)), sub["is_reg"].values, sub["log_total"].values])
y = sub["detected"].values.astype(float)
b, se = logistic(X, y)
z = b/se; p = 2*stats.norm.sf(np.abs(z))
terms = ["intercept", "RegC_HomC", "log_total"]
logit = dict(terms=terms, beta=b.tolist(), se=se.tolist(), z=z.tolist(), p=p.tolist(),
             OR=[math.exp(v) for v in b])
print("\nDepth-adjusted logistic (detected ~ RegC/HomC + log_total):", flush=True)
for t, bb, ss, zz, pp, o in zip(terms, b, se, z, p, logit["OR"]):
    print(f"  {t:10s} beta={bb:+.3f} se={ss:.3f} z={zz:+.2f} p={pp:.4g} OR={o:.3f}", flush=True)

json.dump(dict(fisher=res, per_donor=per_donor, logistic=logit,
               by_state=summary.to_dict(orient="records")),
          open(OUT + r"\GSE324993_replication_results.json", "w"), indent=2)
print("\nDONE", flush=True)
