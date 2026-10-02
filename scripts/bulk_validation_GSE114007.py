# -*- coding: utf-8 -*-
"""Bulk validation of VDR / NLRP3-axis / ligand-receptor findings in an independent
human OA cartilage RNA-seq dataset (GSE114007: 22 OA vs 20 normal)."""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = r"C:\Users\liuha\OpenClawVDR\bulk"
OUT = "C:\\Users\\liuha\\Desktop\\\u5b5f\u5fb7\u5c14\\revision_20261002\\outputs"

oa = pd.read_csv(B + "\\GSE114007_OA_counts.txt.gz", sep="\t", index_col=0)
nm = pd.read_csv(B + "\\GSE114007_normal_counts.txt.gz", sep="\t", index_col=0)
print("OA", oa.shape, "Normal", nm.shape)

data = pd.concat([oa, nm], axis=1)
data = data[~data.index.duplicated()]
X = data.astype(float).replace([np.inf, -np.inf], np.nan).dropna()   # already log2

AXIS = {
 "NLRP3_axis": ["NLRP3","PYCARD","CASP1","IL1B","IL18","GSDMD","TXNIP","AIM2","NLRC4","CASP8"],
 "ligand_receptor": ["SPP1","CD44","THBS1","CD47","FN1","ITGA5","ITGB1","FGF2","FGFR1","APOE","LRP1"],
 "catabolic": ["MMP13","COL10A1","ADAMTS5","IBSP","RUNX2","IL6","CXCL8"],
 "homeostatic": ["ACAN","COL2A1","SOX9","PRG4","COMP","TIMP3"],
}
groups = np.array(["OA"] * oa.shape[1] + ["Normal"] * nm.shape[1])

def zscore_score(genes):
    g = [x for x in genes if x in X.index]
    Z = X.loc[g]
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1, ddof=0).replace(0, 1), axis=0)
    return Z.mean(axis=0).values, g

scores = {}
for k, g in AXIS.items():
    v, used = zscore_score(g)
    scores[k] = v
    print(k, "genes used:", used)

res = []
def cmp(name, vec):
    a = vec[groups == "OA"]; b = vec[groups == "Normal"]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    res.append(dict(feature=name, OA_mean=round(float(a.mean()), 3),
                    Normal_mean=round(float(b.mean()), 3),
                    diff=round(float(a.mean() - b.mean()), 3), p_MWU=float(p)))

cmp("VDR", X.loc["VDR"].values)
for k, v in scores.items():
    cmp(k, v)
for g in ["NLRP3","IL1B","CASP1","PYCARD","GSDMD","SPP1","CD44","THBS1","CD47","MMP13",
          "COL10A1","ACAN","COL2A1","SOX9","CYP24A1"]:
    if g in X.index:
        cmp(g, X.loc[g].values)

resdf = pd.DataFrame(res)
resdf.to_csv(OUT + "\\bulk_GSE114007_OA_vs_normal.csv", index=False)
print("\n" + resdf.to_string(index=False))

cor = []
v = X.loc["VDR"].values
feats = dict(scores)
for g in ["NLRP3","IL1B","CASP1","TXNIP","SPP1","CD44","THBS1","CD47","MMP13","COL10A1",
          "ACAN","COL2A1","SOX9","PRG4","CYP24A1"]:
    if g in X.index:
        feats[g] = X.loc[g].values
for k, arr in feats.items():
    arr = np.asarray(arr)
    for grp in ["OA", "Normal", "All"]:
        m = np.ones(len(v), bool) if grp == "All" else (groups == grp)
        r, p = stats.spearmanr(v[m], arr[m])
        cor.append(dict(feature=k, group=grp, spearman_r=round(float(r), 3), p=float(p), n=int(m.sum())))
cordf = pd.DataFrame(cor)
cordf.to_csv(OUT + "\\bulk_GSE114007_VDR_correlations.csv", index=False)
print("\nVDR correlations (OA only):")
print(cordf[cordf.group == "OA"][["feature","spearman_r","p","n"]].to_string(index=False))

# ---- figure ----
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
a = X.loc["VDR"].values[groups == "OA"]; b = X.loc["VDR"].values[groups == "Normal"]
ax[0].boxplot([a, b], widths=0.6)
ax[0].set_xticks([1, 2]); ax[0].set_xticklabels(["OA", "Normal"])
ax[0].set_ylabel("VDR log2 normalized expression"); ax[0].set_title("A  VDR expression")
u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
ax[0].set_xlabel(f"Mann-Whitney p = {p:.3g}")
a = scores["NLRP3_axis"][groups == "OA"]; b = scores["NLRP3_axis"][groups == "Normal"]
ax[1].boxplot([a, b], widths=0.6)
ax[1].set_xticks([1, 2]); ax[1].set_xticklabels(["OA", "Normal"])
ax[1].set_ylabel("NLRP3-axis score (mean z)"); ax[1].set_title("B  NLRP3-axis activity")
u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
ax[1].set_xlabel(f"Mann-Whitney p = {p:.3g}")
m = groups == "OA"
r, p = stats.spearmanr(X.loc["VDR"].values[m], scores["NLRP3_axis"][m])
ax[2].scatter(X.loc["VDR"].values[m], scores["NLRP3_axis"][m], s=30, color="#2b6cb0")
ax[2].set_xlabel("VDR (log2)"); ax[2].set_ylabel("NLRP3-axis score")
ax[2].set_title(f"C  OA cartilage (n = {int(m.sum())})")
ax[2].text(0.05, 0.93, f"Spearman r = {r:.2f}\np = {p:.3g}", transform=ax[2].transAxes, va="top")
plt.tight_layout()
plt.savefig(OUT + "\\Figure_validation_GSE114007.png", dpi=300)
print("\nfigure saved")
print("DONE_BULK")
