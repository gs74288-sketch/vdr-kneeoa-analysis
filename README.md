# VDR – knee osteoarthritis: analysis code

Code and summary tables supporting the manuscript *"Nominal genetic prioritization of vitamin D
receptor marks a donor-adjusted regulatory chondrocyte state in osteoarthritis"* (revision).

## Data sources (all public)

| Purpose | Source | ID |
|---|---|---|
| Exposure (genetically predicted VDR expression) | eQTLGen whole-blood eQTLs via OpenGWAS | `eqtl-a-ENSG00000111424` |
| Outcome (knee OA) | Tachmazidou et al. 2019 via OpenGWAS / GWAS Catalog | `ebi-a-GCST007090` / GCST007090 (24,955 cases / 378,169 controls) |
| Single-cell atlas (primary) | GEO | GSE104782 |
| Single-cell atlas (independent replication) | GEO | GSE152805 (OA knee cartilage, 6 samples, 3 donors) |

## Analysis steps

1. **Target-focused MR screen + BH-FDR** — `coloc_steiger_final.R`, `cis_mr_and_recoloc.R`
2. **Colocalization** (`coloc::coloc.abf`) over chr12:46.8–48.9 Mb (GRCh37), 6,988 harmonised SNPs — `coloc_steiger_final.R`
3. **Steiger directionality** (binary outcome corrected with `get_r_from_lor`) — `steiger_only.R`
4. **Cis-restricted MR sensitivity** (distance clumping; OpenGWAS LD-clumping endpoint was unavailable) — `cis_mr_local.R`
5. **Independent scRNA-seq replication** — `replicate_GSE152805.R` + `gse152805_analysis.R`
6. OpenGWAS API helpers — `fetch_region.py`, `fetch_opengwas.py`

R 4.4.3; TwoSampleMR, coloc, Seurat 5, ieugwasr, data.table.

## Key results

- **Colocalization** (VDR eQTL vs knee OA): PP.H3 = 0.971, PP.H4 = 1.9 × 10⁻⁴ (no shared causal variant).
- **Steiger**: correct causal direction = TRUE; P = 2.7 × 10⁻¹⁴⁵ (snp_r2 exposure 0.025 vs outcome 1.1 × 10⁻⁴).
- **Cis-restricted MR** (distance-clumped): IVW OR = 0.876 (100 kb), OR = 0.859 (10 kb); top cis SNP
  rs7975232 alone not significant (OR 0.960, P = 0.41). LD clumping not applied — estimates anticonservative.
- **Independent replication (GSE152805)**: VDR detection did **not** replicate enrichment in homeostatic
  chondrocytes (HTC 6.9% > FC 5.3% > HomC 4.9% > EC 4.4%; Fisher OR 1.01, P = 0.94; donor/depth-adjusted
  OR 0.91, P = 0.20).

`tables/` contains the summary outputs; `scripts/` the code.
