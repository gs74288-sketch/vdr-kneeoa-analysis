# -*- coding: utf-8 -*-
"""Project the knee OA GWAS onto chondrocyte states with scDRS.

Input : data/GSE324993.h5ad        (AnnData; obs: sample, group, cell_type)
        tables/GCST007090_gene_z.tsv
Output: tables/scDRS_cell_scores.csv
        tables/scDRS_by_state.csv
        tables/scDRS_per_donor.csv

Notes
-----
* scDRS normalises scores across cells, so absolute scores are not comparable
  between datasets; differences are interpreted within this atlas only.
* With three donors per group a case-control comparison is not identifiable;
  the informative test is the within-donor contrast between RegC and other
  states, aggregated across donors.
"""
import os

import numpy as np
import pandas as pd
import anndata as ad
import scdrs
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "tables")

H5 = os.path.join(DATA, "GSE324993.h5ad")
Z = os.path.join(OUT, "GCST007090_gene_z.tsv")

z = pd.read_csv(Z, sep="\t").dropna()
z = z[~z["GENE"].duplicated()]
gene_list = z["GENE"].values
gene_weight = z["Z"].values.astype(float)
print("gene set:", len(gene_list))

adata = ad.read_h5ad(H5)
print("loaded:", adata.shape)
scdrs.preprocess(adata, n_mean_bin=20, n_var_bin=20, copy=False)

df_score = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000,
                            return_ctrl_raw_score=False, return_ctrl_norm_score=False)
df_score.to_csv(os.path.join(OUT, "scDRS_cell_scores.csv"))

obs = adata.obs.join(df_score, how="inner")
col = "norm_score"

# library-size residual (scores are not fully depth-independent)
X = adata.X.tocsr()
obs["log_counts"] = np.log1p(np.asarray(X.sum(axis=1)).ravel())
b1, b0 = np.polyfit(obs["log_counts"], obs[col], 1)
obs["resid"] = obs[col] - (b1 * obs["log_counts"] + b0)
obs["is_reg"] = obs["cell_type"].astype(str).str.startswith("RegC")

by_state = (obs.assign(state=obs["cell_type"].astype(str))
               .groupby("state")[col]
               .agg(n="size", mean="mean", median="median")
               .reset_index()
               .sort_values("mean", ascending=False))
by_state.to_csv(os.path.join(OUT, "scDRS_by_state.csv"), index=False)
print(by_state.to_string(index=False))

rows = []
for s, sub in obs.groupby("sample"):
    a = sub.loc[sub.is_reg, "resid"]
    b = sub.loc[~sub.is_reg, "resid"]
    if len(a) > 20 and len(b) > 20:
        rows.append(dict(donor=s, group=sub["group"].iloc[0],
                         n_reg=int(len(a)), n_other=int(len(b)),
                         delta=float(a.mean() - b.mean())))
per = pd.DataFrame(rows)
per.to_csv(os.path.join(OUT, "scDRS_per_donor.csv"), index=False)
print(per.to_string(index=False))

d = per["delta"].values
t, p = stats.ttest_1samp(d, 0.0)
print("\ndonor-level: mean %+.4f, SD %.4f, t = %.2f, P = %.3f, positive %d/%d" % (
    d.mean(), d.std(ddof=1), t, p, int((d > 0).sum()), len(d)))
print("DONE")
