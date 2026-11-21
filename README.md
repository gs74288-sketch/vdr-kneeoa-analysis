# VDR and knee osteoarthritis: analysis code

Code and summary tables supporting the manuscript *"The knee osteoarthritis signal at 12q13 maps to COL2A1, not the vitamin D receptor"*.

## Data sources (all public)

| Purpose | Source | ID |
|---|---|---|
| Exposure (genetically predicted VDR expression) | eQTLGen whole-blood eQTLs via OpenGWAS | `eqtl-a-ENSG00000111424` |
| Outcome (knee OA, primary) | 2025 multi-ancestry knee OA GWAS via GWAS Catalog | `GCST90566800` (161,228 European cases / 965,985 European controls) |
| Outcome (knee OA, exploratory screen and FinnGen comparison) | Tachmazidou et al. 2019 via OpenGWAS / GWAS Catalog | `ebi-a-GCST007090` / GCST007090 (24,955 cases / 378,169 controls) |
| Single-cell atlas (primary) | GEO | GSE104782 |
| Single-cell atlas (replication) | GEO | GSE152805 (OA knee cartilage, 6 samples, 3 donors) |
| Single-cell atlas (replication) | GEO | GSE324993 (OA and healthy cartilage, 44,348 cells, 6 donors) |
| Bulk cartilage | GEO | GSE114007 |

## Analysis steps

1. **Exploratory target-focused MR screen + BH-FDR** — `coloc_steiger_final.R`, `cis_mr_and_recoloc.R`
2. **Gene attribution in the 12q13 interval** (gene bodies vs. genome-wide significant variants) — `_gene_attrib` logic in `coloc_GCST90566800.R` companion script
3. **Colocalisation** (`coloc::coloc.abf`) over chr12:46.8-48.9 Mb (GRCh38), 6,496 harmonised SNPs against GCST90566800 — `coloc_GCST90566800.R`
4. **Colocalisation window sensitivity** (±10 to ±100 kb around the expression and disease peaks; p12 1e-6 to 1e-4) — `coloc_windows.R`
5. **Steiger directionality** (binary outcome corrected with `get_r_from_lor`) — `steiger_only.R`
6. **Cis-restricted MR, LD-clumped** against the 1000 Genomes European reference panel (`r2 < 0.001`) using LDlink/LDlinkR — `cis_mr_local.R`
7. **Sensitivity: genome-wide significant variants only** (P < 5e-8), fixed and random effects — `cis_mr_GWonly.R`
8. **Cis-restricted sensitivity with distance clumping** — `cis_mr_local.R`
9. **Independent scRNA-seq replication, GSE152805** — `replicate_GSE152805.R` + `gse152805_analysis.R`
10. **Independent scRNA-seq replication, GSE324993** (donor-stratified; VDR detection by state) — `replicate_GSE324993.py`
11. **Gene-level z-scores for scDRS** (GTF gene bodies +/- 10 kb, min p per gene) — `build_gene_z.py`
12. **Single-cell disease relevance scoring** — `run_scdrs.py`
13. **Bulk validation, GSE114007** — `bulk_validation_GSE114007.py`
14. OpenGWAS API helpers — `fetch_region.py`, `fetch_opengwas.py`

R 4.4.3; TwoSampleMR, coloc, Seurat 5, ieugwasr, data.table, scTenifoldKnk. Python 3.11; scdrs 1.0.2, scanpy, anndata.

## Key results

- **The 12q13 knee OA signal belongs to COL2A1, not VDR.** Against GCST90566800 (161,228 European cases), 80 variants in chr12:46.8-48.9 Mb reached genome-wide significance and every one lies closer to COL2A1 (47,972,967-48,004,554) than to VDR (47,841,537-47,943,048). None of the 814 variants within the VDR gene body was significant (best P = 3.3 × 10^-3); 20 within COL2A1 were (best P = 2.0 × 10^-11). The VDR eQTL peak rs7975232 (47,845,054; P = 6.5 × 10^-56) lies 165 kb from the disease peak rs11168351 (48,009,982; P = 1.5 × 10^-12) and has no OA association (P = 0.58).
- **Colocalisation**: PP.H3 = 1.00, PP.H4 = 2.2 × 10^-6 over the full interval, unchanged for p12 from 1 × 10^-6 to 1 × 10^-4. In ±10 kb and ±100 kb windows around the expression peak, where knee OA has no signal at all, the posterior returns association with expression alone (PP.H1 = 0.992 and 0.980).
- **Cis-restricted MR is not interpretable at this locus.** Four LD-clumped instruments were defined from expression data alone; three were present in GCST90566800 (rs191597566 absent). rs4506714, which sits between TMEM106C and COL2A1 and is itself associated with knee OA (P = 1.1 × 10^-7), gives Wald OR 0.876 (0.835-0.920), whereas the strongest expression variant rs7975232 gives 1.013. Random-effects IVW OR = 0.954 (0.858-1.061), P = 0.38; fixed-effect 0.954 (0.923-0.985), P = 0.004; **Cochran Q = 21.5 on 2 df (P = 2.2 × 10^-5)**. The heterogeneity shows the instruments are not estimating a single quantity, so no causal interpretation is offered.
- **Steiger directionality**: correct causal direction = TRUE; P = 2.7 × 10^-145.
- **FinnGen consistency check** (FinnGen also contributes to GCST90566800, so this is not an independent replication): exploratory instrument set IVW OR = 0.965 (P = 0.34); knee-surgery endpoint OR = 0.933 (P = 0.11).
- **VDR detection across three atlases** (donor- and depth-adjusted, regulatory vs other states): GSE104782 OR = 1.89 (1.32-2.69); GSE152805 OR = 1.50 (1.09-2.06); GSE324993 OR = 1.13 (1.07-1.19), positive in all six donors. State definitions differ between atlases, so no pooled estimate is treated as primary; an exploratory random-effects meta-analysis gives OR = 1.42 (1.02-1.97) with I^2 = 81% and a 95% prediction interval of 0.38-5.33.
- **Disease relevance (scDRS)** on GSE324993: highest in RepC (mean normalised score +0.49); the merged regulatory compartment (16,905 cells) scored +0.06 versus -0.04 in other states; within-donor difference positive in 5 of 6 donors (mean +0.17, one-sample t = 2.77, P = 0.040). A case-control comparison was not interpretable with three donors per group and is not claimed.
- **Bulk validation (GSE114007)**: VDR upregulated in OA cartilage (log2 3.93 vs 3.58, P = 8.5 × 10^-3) and correlated with CYP24A1 (Spearman r = 0.50, P = 0.017).

## Pitfalls encountered (worth knowing before re-running)

- **Truncated downloads fail silently.** `GSE324993_counts_matrix.csv.gz` was downloaded at 22.4 MB instead of 113.4 MB; `zlib.decompressobj` stopped at the truncation point without raising, and the VDR row was never read, producing all-zero detection. `replicate_GSE324993.py` now checks `decompressor.eof` and aborts. Any streaming parser should do the same.
- **Genome build mismatch.** EBI GWAS Catalog harmonised files are GRCh38. Using a GRCh37 GTF shifted gene assignments by up to ~1.5 Mb and put GDF5 at rank 38,475 (true Z ≈ +9.9, rank 1). Check a known locus before trusting gene-level output.
- **Two analyses that did not work here** and are not reported in the manuscript: bulk deconvolution (NNLS and MuSiC both fail because GSE114007 is length-normalised bulk while the reference is 10x UMI; MuSiC reconstruction R^2 was 0.03-0.09), and scTenifoldKnk virtual knockout (control knockouts of unrelated genes produced equal or larger numbers of significant genes, so the ranking is not target-specific).
- **Gene-dense regions break cis-eQTL MR.** Four LD-clumped VDR instruments were defined from expression alone; one of them (rs4506714) tags the adjacent COL2A1 signal and carries a genome-wide significant disease association. Against the smaller GCST007090 the instrument set looked homogeneous (Q = 5.04 on 2 df, P = 0.080); against GCST90566800 the same set is grossly heterogeneous (Q = 21.5 on 2 df, P = 2.2 × 10^-5). Weak outcome data cannot resolve instrument heterogeneity, which is why a null-looking but homogeneous IVW is not evidence of absence.
- **Per-variant sample size varies within one GWAS file.** GCST90566800 reports several rows per position with different case counts; we kept the row with the largest case count per rsID and derived standard errors from the published 95% confidence intervals because the `standard_error` column was empty.
- **Whole-file download is avoidable.** GCST90566800 is 1.12 GB over a slow link; reading the `.tbi` index and fetching only the chr12 chunks took 951 KB. See the tabix/BGZF reader used to build `tables/coloc_harmonised_input_GCST90566800.csv`.

`tables/` contains the summary outputs; `scripts/` the code.
