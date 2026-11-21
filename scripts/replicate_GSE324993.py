# -*- coding: utf-8 -*-
"""Directional replication of VDR detection enrichment in regulatory chondrocyte
states, independent atlas GSE324993 (OA + healthy cartilage scRNA-seq).

Handles a gzip stream whose final trailer may be missing at source (zlib
streaming) and VERIFIES that the stream completed: a truncated download makes
zlib stop early without error, which previously produced all-zero detection.

Input : data/GSE324993_counts_matrix.csv.gz   (genes x cells, CSV inside gzip)
        data/GSE324993_cell_annotation.txt.gz
Output: tables/GSE324993_VDR_by_state.csv
        tables/GSE324993_VDR_per_donor.csv
        tables/GSE324993_replication.json
"""
import gzip
import json
import math
import os
import zlib

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "tables")
os.makedirs(OUT, exist_ok=True)

MAT = os.path.join(DATA, "GSE324993_counts_matrix.csv.gz")
ANN = os.path.join(DATA, "GSE324993_cell_annotation.txt.gz")

# ---- cell annotation ----
raw2type = {}
with gzip.open(ANN, "rt", encoding="utf-8", errors="replace") as f:
    started = False
    for line in f:
        line = line.rstrip("\n")
        if not started:
            if line.startswith("Cell_Label"):
                started = True
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            raw2type[parts[0]] = parts[1]

CANON = {"regc": "RegC", "homc": "HomC", "htc": "HTC", "prehtc": "preHTC", "proc": "ProC",
         "ec": "EC", "fc": "FC", "mtc": "MTC", "cpc": "CPC", "prefc": "PreFC", "repc": "RepC"}


def canon(ct):
    s = ct.split("(")[0].strip().lower()
    for k, v in CANON.items():
        if s.startswith(k):
            return v
    return ct


# ---- stream matrix, checking completeness ----
GENES = ["VDR"]
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
        buf += dec.decompress(chunk).decode("utf-8", errors="replace")
        *lines, buf = buf.split("\n")
        for line in lines:
            line = line.rstrip("\r")
            if not line:
                continue
            if cells is None:
                cells = line.split(",")
                totals = np.zeros(len(cells), dtype=np.int64)
                continue
            i = line.find(",")
            if i < 0:
                continue
            g = line[:i]
            vals = np.fromiter(map(int, line[i + 1:].split(",")),
                               dtype=np.int64, count=len(cells))
            if len(vals) != len(cells):
                continue
            totals += vals
            if g in GENES:
                gene_rows[g] = vals.astype(np.float32)
            nrow += 1

# the critical check that was missing before
if not dec.eof:
    raise RuntimeError(
        "gzip stream ended prematurely (decompressor.eof is False). "
        "The count matrix is truncated; re-download the file before analysing. "
        "Genes parsed so far: %d" % nrow)

print("cells:", len(cells), "| genes parsed:", nrow, "| VDR captured:", "VDR" in gene_rows)

vdr = gene_rows.get("VDR")
if vdr is None:
    raise RuntimeError("VDR not found in the count matrix")

df = pd.DataFrame({"cell": cells, "vdr": vdr, "total": totals})
df["sample"] = df["cell"].map(lambda c: c.split("_")[0])
df["group"] = np.where(df["sample"].str.startswith("OA"), "OA",
                       np.where(df["sample"].str.startswith("HC"), "HC", "NA"))
df["state"] = df["cell"].map(lambda c: canon(raw2type.get(c, "NA")))
df["detected"] = (df["vdr"] > 0).astype(int)

summary = (df.groupby("state")
             .agg(n=("cell", "size"), detected=("detected", "sum"),
                  mean_total=("total", "mean"))
             .reset_index())
summary["frac"] = summary["detected"] / summary["n"]
summary = summary.sort_values("frac", ascending=False)
summary.to_csv(os.path.join(OUT, "GSE324993_VDR_by_state.csv"), index=False)
print(summary.to_string(index=False))

# ---- per donor ----
per = []
for s in sorted(df["sample"].unique()):
    sub = df[df["sample"] == s]
    reg = sub[sub["state"] == "RegC"]
    oth = sub[sub["state"] != "RegC"]
    if len(reg) and len(oth):
        per.append(dict(sample=s, group=sub["group"].iloc[0],
                        reg_n=int(len(reg)), reg_frac=float(reg["detected"].mean()),
                        other_n=int(len(oth)), other_frac=float(oth["detected"].mean()),
                        delta=float(reg["detected"].mean() - oth["detected"].mean())))
pd.DataFrame(per).to_csv(os.path.join(OUT, "GSE324993_VDR_per_donor.csv"), index=False)
for r in per:
    print(r)
print("donors with positive delta: %d/%d" % (sum(1 for r in per if r["delta"] > 0), len(per)))

# ---- donor- and depth-adjusted logistic ----
states_ok = df["state"] != "NA"
sub = df[states_ok].copy()
sub["is_reg"] = sub["state"].eq("RegC").astype(float)
sub["log_total"] = np.log1p(sub["total"])
donors = sorted(sub["sample"].unique())
X = [np.ones(len(sub)), sub["is_reg"].values, sub["log_total"].values]
names = ["intercept", "RegC", "log_total"]
for d in donors[:-1]:
    X.append((sub["sample"] == d).astype(float).values)
    names.append("donor_" + d)
X = np.column_stack(X)
y = sub["detected"].values.astype(float)

beta = np.zeros(X.shape[1])
for _ in range(80):
    eta = X @ beta
    mu = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
    w = np.clip(mu * (1 - mu), 1e-9, None)
    z = eta + (y - mu) / w
    XtW = X.T * w
    bn = np.linalg.solve(XtW @ X, XtW @ z)
    if np.max(np.abs(bn - beta)) < 1e-10:
        beta = bn
        break
    beta = bn
eta = X @ beta
mu = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
w = np.clip(mu * (1 - mu), 1e-9, None)
se = np.sqrt(np.diag(np.linalg.pinv((X * w[:, None]).T @ X)))
i = names.index("RegC")
print("\nRegC (donor + depth adjusted): beta=%.4f se=%.4f OR=%.3f (%.3f-%.3f) P=%.3g" % (
    beta[i], se[i], math.exp(beta[i]),
    math.exp(beta[i] - 1.96 * se[i]), math.exp(beta[i] + 1.96 * se[i]),
    2 * stats.norm.sf(abs(beta[i] / se[i]))))

json.dump(dict(by_state=summary.to_dict(orient="records"), per_donor=per,
               adjusted=dict(term="RegC", beta=float(beta[i]), se=float(se[i]),
                             OR=float(math.exp(beta[i])),
                             ci_low=float(math.exp(beta[i] - 1.96 * se[i])),
                             ci_high=float(math.exp(beta[i] + 1.96 * se[i])))),
          open(os.path.join(OUT, "GSE324993_replication.json"), "w"), indent=2)
print("\nDONE")
