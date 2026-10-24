# VDR and knee osteoarthritis: analysis code

Code and summary tables supporting the manuscript *"The vitamin D receptor in knee osteoarthritis: confounded association, no causal genetic support, and reactive upregulation in cartilage"*.

## Data sources (all public)

| Purpose | Source | ID |
|---|---|---|
| Exposure (genetically predicted VDR expression) | eQTLGen whole-blood eQTLs via OpenGWAS | `eqtl-a-ENSG00000111424` |
| Outcome (knee OA) | Tachmazidou et al. 2019 via OpenGWAS / GWAS Catalog | `ebi-a-GCST007090` / GCST007090 (24,955 cases / 378,169 controls) |
| Single-cell atlas (primary) | GEO | GSE104782 |
| Single-cell atlas (independent replication) | GEO | GSE152805 (OA knee cartilage, 6 samples, 3 donors) |

## Analysis steps

1. **Exploratory target-focused MR screen + BH-FDR** — `coloc_steiger_final.R`, `cis_mr_and_recoloc.R`
2. **Colocalisation** (`coloc::coloc.abf`) over chr12:46.8-48.9 Mb (GRCh37), 6,988 harmonised SNPs — `coloc_steiger_final.R`
3. **Steiger directionality** (binary outcome corrected with `get_r_from_lor`) — `steiger_only.R`
4. **Primary cis-restricted MR, LD-clumped** against the 1000 Genomes European reference panel (`r2 < 0.001`) using LDlink/LDlinkR — `cis_mr_local.R`
5. **Cis-restricted sensitivity with distance clumping** (the OpenGWAS LD-clumping endpoint was unavailable at the time) — `cis_mr_local.R`
6. **Independent scRNA-seq replication** — `replicate_GSE152805.R` + `gse152805_analysis.R`
7. OpenGWAS API helpers — `fetch_region.py`, `fetch_opengwas.py`

R 4.4.3; TwoSampleMR, coloc, Seurat 5, ieugwasr, data.table.

## Key results

- **Primary cis-restricted MR** (LD-clumped, 4 instruments; 3 reach genome-wide significance): IVW OR = 0.915 (95% CI 0.831-1.007), P = 0.069; weighted median OR = 0.913, P = 0.036; no heterogeneity (Cochran Q P = 0.12) or directional pleiotropy (intercept P = 0.37). The single strongest cis variant, rs7975232, analysed alone: OR = 0.960 (0.870-1.059), P = 0.414.
- **Cis-restricted distance-clumped sensitivity**: IVW OR = 0.876 (100 kb) and OR = 0.859 (10 kb). Reported for transparency; distance clumping does not account for linkage disequilibrium between retained variants, so these estimates are anticonservative and are not the primary result.
- **Colocalisation** (VDR eQTL vs knee OA): PP.H3 = 0.971, PP.H4 = 1.9 × 10^-4 (no shared causal variant).
- **Steiger**: correct causal direction = TRUE; P = 2.7 × 10^-145 (snp_r2 exposure 0.025 vs outcome 1.1 × 10^-4).
- **FinnGen cis replication**: OR = 0.984 (P = 0.55); knee-surgery endpoint OR = 0.968 (P = 0.38).
- **Independent replication (GSE152805)**: VDR detection did **not** replicate enrichment in homeostatic chondrocytes (HTC 6.9% > FC 5.3% > HomC 4.9% > EC 4.4%; Fisher OR 1.01, P = 0.94; donor/depth-adjusted OR 0.91, P = 0.20).

`tables/` contains the summary outputs; `scripts/` the code.
