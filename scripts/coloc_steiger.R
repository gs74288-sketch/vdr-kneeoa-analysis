# Colocalization + Steiger directionality for VDR -> knee OA
# Inputs produced by filter_gwas_kneeoa.py and filter_eqtlgen_vdr.py
suppressMessages({
  library(data.table)
  library(TwoSampleMR)
  library(coloc)
})

base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
gw_f <- file.path(base, "data/GCST007090_chr12_VDR_region.tsv")
eq_f <- file.path(base, "data/eQTLGen_VDR_cis.tsv")
outdir <- file.path(base, "outputs")
dir.create(outdir, showWarnings = FALSE)

N_EQTL <- 31684      # eQTLGen
NCASE_OA <- 24955    # Tachmazidou 2019 knee OA cases
NCONTROL_OA <- 378169
N_OA <- NCASE_OA + NCONTROL_OA
S_OA <- NCASE_OA / N_OA

gw <- fread(gw_f)
eq <- fread(eq_f)
cat("GWAS region rows:", nrow(gw), " eQTLGen VDR rows:", nrow(eq), "\n")
cat("GWAS cols:", paste(names(gw), collapse=","), "\n")
cat("eQTL cols:", paste(names(eq), collapse=","), "\n")

# --- harmonise by rsID ---
gw2 <- gw[!is.na(hm_rsid) & hm_rsid != "" & hm_rsid != "."]
gw2 <- gw2[, .(snp = hm_rsid,
               a1 = toupper(hm_effect_allele),
               a2 = toupper(hm_other_allele),
               beta_gwas = as.numeric(hm_beta),
               se_gwas = as.numeric(standard_error),
               p_gwas = as.numeric(p_value),
               eaf_gwas = as.numeric(hm_effect_allele_frequency),
               pos = as.integer(hm_pos))]
gw2 <- gw2[!is.na(beta_gwas) & !is.na(se_gwas) & !is.na(snp)]
gw2 <- gw2[!duplicated(snp)]

eq2 <- as.data.table(eq)
eq2 <- eq2[, .(snp = SNP,
               ea = toupper(AssessedAllele),
               oa = toupper(OtherAllele),
               z = as.numeric(Zscore),
               p_eqtl = as.numeric(Pvalue),
               n_eqtl = as.integer(NrSamples))]
eq2 <- eq2[!is.na(z) & !is.na(snp)]
eq2 <- eq2[!duplicated(snp)]

m <- merge(eq2, gw2, by = "snp")
cat("merged SNPs:", nrow(m), "\n")

# align eQTL effect allele to GWAS effect allele (a1)
flip <- m$ea != m$a1
# if alleles are swapped (ea==a2), flip z sign
swap <- m$ea == m$a2
m$z_adj <- ifelse(swap, -m$z, m$z)
# MAF from GWAS eaf (aligned to GWAS effect allele a1)
m$maf <- pmin(m$eaf_gwas, 1 - m$eaf_gwas)
m <- m[!is.na(m$maf) & m$maf > 0 & m$maf < 1]

n_e <- ifelse(is.na(m$n_eqtl), N_EQTL, m$n_eqtl)
denom <- sqrt(2 * m$maf * (1 - m$maf) * (n_e + m$z_adj^2))
m$beta_eqtl <- m$z_adj / denom
m$se_eqtl <- 1 / denom

# --- COLOC (coloc.abf, single causal variant) ---
d1 <- list(snp = m$snp, beta = m$beta_eqtl, varbeta = m$se_eqtl^2,
           type = "quant", N = N_EQTL, MAF = m$maf, sdY = 1)
d2 <- list(snp = m$snp, beta = m$beta_gwas, varbeta = m$se_gwas^2,
           type = "cc", N = N_OA, s = S_OA, MAF = m$maf)
res <- coloc.abf(d1, d2)
summary_df <- as.data.frame(t(res$summary))
fwrite(summary_df, file.path(outdir, "coloc_VDR_kneeOA_summary.csv"))
pp <- res$results
pp <- pp[order(-pp$SNP.PP.H4), ]
fwrite(pp, file.path(outdir, "coloc_VDR_kneeOA_perSNP.csv"))
cat("COLOC summary:\n"); print(res$summary)

# --- Steiger directionality (TwoSampleMR) ---
instr <- fread(file.path(base, "..", "SuppTable3_VDR_SNPs_Details_and_F_stats.csv"))
instr <- as.data.table(instr)
setnames(instr, c("SNP","chr.exposure","pos.exposure","effect_allele.exposure",
                  "other_allele.exposure","eaf.exposure","beta.exposure",
                  "se.exposure","pval.exposure","F_statistic"))
expo <- data.frame(
  SNP = instr$SNP,
  beta.exposure = as.numeric(instr$beta.exposure),
  se.exposure = as.numeric(instr$se.exposure),
  effect_allele.exposure = instr$effect_allele.exposure,
  other_allele.exposure = instr$other_allele.exposure,
  eaf.exposure = as.numeric(instr$eaf.exposure),
  samplesize.exposure = N_EQTL,
  units.exposure = "SD",
  id.exposure = "VDR_eqtlgen"
)
out <- data.frame(
  SNP = m$snp,
  beta.outcome = m$beta_gwas,
  se.outcome = m$se_gwas,
  effect_allele.outcome = m$a1,
  other_allele.outcome = m$a2,
  eaf.outcome = m$eaf_gwas,
  samplesize.outcome = N_OA,
  ncase.outcome = NCASE_OA,
  ncontrol.outcome = NCONTROL_OA,
  units.outcome = "log odds",
  id.outcome = "knee_OA"
)
dat <- merge(expo, out, by = "SNP")
dat <- harmonise_data(expo, out, action = 2)
cat("harmonised instruments:", nrow(dat), " kept:", sum(dat$mr_keep), "\n")
st <- tryCatch(directionality_test(dat), error = function(e) { cat("Steiger error:", conditionMessage(e), "\n"); NULL })
if (!is.null(st)) {
  fwrite(as.data.table(st), file.path(outdir, "steiger_VDR_kneeOA.csv"))
  print(st)
}
cat("DONE coloc+steiger\n")
