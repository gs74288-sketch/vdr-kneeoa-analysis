# VDR and knee osteoarthritis: analysis code

Code and summary tables supporting the manuscript *"No evidence of a shared causal variant for vitamin D receptor expression and knee osteoarthritis: a triangulation of colocalisation, population and cartilage evidence"*.

## Data sources (all public)

| Purpose | Source | ID |
|---|---|---|
| Exposure (genetically predicted VDR expression) | eQTLGen whole-blood eQTLs via OpenGWAS | `eqtl-a-ENSG00000111424` |
| Outcome (knee OA) | Tachmazidou et al. 2019 via OpenGWAS / GWAS Catalog | `ebi-a-GCST007090` / GCST007090 (24,955 cases / 378,169 controls) |
| Single-cell atlas (primary) | GEO | GSE104782 |
| Single-cell atlas (replication) | GEO | GSE152805 (OA knee cartilage, 6 samples, 3 donors) |
| Single-cell atlas (replication) | GEO | GSE324993 (OA and healthy cartilage, 44,348 cells, 6 donors) |
| Bulk cartilage | GEO | GSE114007 |

## Analysis steps

1. **Exploratory target-focused MR screen + BH-FDR** — `coloc_steiger_final.R`, `cis_mr_and_recoloc.R`
2. **Colocalisation** (`coloc::coloc.abf`) over chr12:46.8-48.9 Mb, 6,988 harmonised SNPs — `coloc_steiger_final.R`
3. **Steiger directionality** (binary outcome corrected with `get_r_from_lor`) — `steiger_only.R`
4. **Primary cis-restricted MR, LD-clumped** against the 1000 Genomes European reference panel (`r2 < 0.001`) using LDlink/LDlinkR — `cis_mr_local.R`
5. **Sensitivity: genome-wide significant variants only** (P < 5e-8), fixed and random effects — `cis_mr_GWonly.R`
6. **Cis-restricted sensitivity with distance clumping** — `cis_mr_local.R`
7. **Independent scRNA-seq replication, GSE152805** — `replicate_GSE152805.R` + `gse152805_analysis.R`
8. **Independent scRNA-seq replication, GSE324993** (donor-stratified; VDR detection by state) — `replicate_GSE324993.py`
9. **Gene-level z-scores for scDRS** (GTF gene bodies +/- 10 kb, min p per gene) — `build_gene_z.py`
10. **Single-cell disease relevance scoring** — `run_scdrs.py`
11. **Bulk validation, GSE114007** — `bulk_validation_GSE114007.py`
12. OpenGWAS API helpers — `fetch_region.py`, `fetch_opengwas.py`

R 4.4.3; TwoSampleMR, coloc, Seurat 5, ieugwasr, data.table, scTenifoldKnk. Python 3.11; scdrs 1.0.2, scanpy, anndata.

## Key results

- **Colocalisation** (VDR eQTL vs knee OA): PP.H3 = 0.971, PP.H4 = 1.9 × 10^-4 — distinct causal variants, so the expression signal is not carried by a disease variant.
- **Primary cis-restricted MR** (LD-clumped, 4 instruments; 3 reach genome-wide significance): IVW OR = 0.915 (95% CI 0.831-1.007), P = 0.069. The strongest instrument, rs7975232, is itself associated with bone mineral density, so the exclusion restriction cannot be assumed and the estimate is supporting rather than decisive.
- **Sensitivity: genome-wide significant variants only** (3 of 4). Point estimate unchanged (OR = 0.918), interval depends on the heterogeneity model: random effects 0.822-1.026 (P = 0.13), fixed effect 0.856-0.984 (P = 0.016). Cochran Q = 5.04 on 2 df (P = 0.080).
- **Cis-restricted distance-clumped sensitivity**: IVW OR = 0.876 (100 kb) and 0.859 (10 kb); anticonservative, not the primary result.
- **Steiger**: correct causal direction = TRUE; P = 2.7 × 10^-145.
- **FinnGen cis replication**: OR = 0.984 (P = 0.55); knee-surgery endpoint OR = 0.968 (P = 0.38).
- **VDR detection across three atlases** (donor- and depth-adjusted, regulatory vs other states): GSE104782 OR = 1.89 (1.32-2.69); GSE152805 OR = 1.50 (1.09-2.06); GSE324993 OR = 1.13 (1.07-1.19), positive in all six donors. State definitions differ between atlases, so no pooled estimate is treated as primary; an exploratory random-effects meta-analysis gives OR = 1.42 (1.02-1.97) with I^2 = 81% and a 95% prediction interval of 0.38-5.33.
- **Disease relevance (scDRS)** on GSE324993: highest in RepC (mean normalised score +0.49); the merged regulatory compartment (16,905 cells) scored +0.06 versus -0.04 in other states; within-donor difference positive in 5 of 6 donors (mean +0.17, one-sample t = 2.77, P = 0.040). A case-control comparison was not interpretable with three donors per group and is not claimed.
- **Bulk validation (GSE114007)**: VDR upregulated in OA cartilage (log2 3.93 vs 3.58, P = 8.5 × 10^-3) and correlated with CYP24A1 (Spearman r = 0.50, P = 0.017).

## Pitfalls encountered (worth knowing before re-running)

- **Truncated downloads fail silently.** `GSE324993_counts_matrix.csv.gz` was downloaded at 22.4 MB instead of 113.4 MB; `zlib.decompressobj` stopped at the truncation point without raising, and the VDR row was never read, producing all-zero detection. `replicate_GSE324993.py` now checks `decompressor.eof` and aborts. Any streaming parser should do the same.
- **Genome build mismatch.** EBI GWAS Catalog harmonised files are GRCh38. Using a GRCh37 GTF shifted gene assignments by up to ~1.5 Mb and put GDF5 at rank 38,475 (true Z ≈ +9.9, rank 1). Check a known locus before trusting gene-level output.
- **Two analyses that did not work here** and are not reported in the manuscript: bulk deconvolution (NNLS and MuSiC both fail because GSE114007 is length-normalised bulk while the reference is 10x UMI; MuSiC reconstruction R^2 was 0.03-0.09), and scTenifoldKnk virtual knockout (control knockouts of unrelated genes produced equal or larger numbers of significant genes, so the ranking is not target-specific).

`tables/` contains the summary outputs; `scripts/` the code.
