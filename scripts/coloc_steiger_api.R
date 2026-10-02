# Colocalization + Steiger directionality for VDR -> knee OA, via OpenGWAS API
suppressMessages({
  library(ieugwasr); library(TwoSampleMR); library(data.table); library(coloc)
})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
outdir <- file.path(base, "outputs"); dir.create(outdir, showWarnings = FALSE)
N_EQTL <- 31684; NCASE <- 24955; NCONTROL <- 378169; N_OA <- NCASE + NCONTROL

cat("== fetching VDR cis-eQTL (eqtl-a-ENSG00000111424) ==\n")
eq <- tryCatch(as.data.table(tophits(id = "eqtl-a-ENSG00000111424", pval = 1)),
               error = function(e){cat("tophits error:", conditionMessage(e), "\n"); NULL})
if (is.null(eq)) quit(status = 1)
cat("rows:", nrow(eq), "\ncols:", paste(names(eq), collapse=","), "\n")
print(head(eq, 3))
fwrite(eq, file.path(outdir, "VDR_cis_eQTL_opengwas.csv"))

snps <- unique(eq$rsid)
cat("== fetching knee OA outcome for", length(snps), "SNPs ==\n")
oa <- tryCatch(as.data.table(associations(variants = snps, id = "ebi-a-GCST007090")),
               error = function(e){cat("assoc error:", conditionMessage(e), "\n"); NULL})
if (is.null(oa)) quit(status = 1)
cat("outcome rows:", nrow(oa), "\ncols:", paste(names(oa), collapse=","), "\n")
fwrite(oa, file.path(outdir, "kneeOA_region_opengwas.csv"))

# ---- harmonise by allele ----
e <- eq[, .(snp=rsid, ea1=toupper(ea), oa1=toupper(oa), eaf1=as.numeric(eaf),
            b1=as.numeric(beta), se1=as.numeric(se), p1=as.numeric(p))]
o <- oa[, .(snp=rsid, ea2=toupper(ea), oa2=toupper(oa), eaf2=as.numeric(eaf),
            b2=as.numeric(beta), se2=as.numeric(se), p2=as.numeric(p))]
m <- merge(e, o, by="snp")
m <- m[ea1 != oa1]
same <- m$ea1 == m$ea2 & m$oa1 == m$oa2
swap <- m$ea1 == m$oa2 & m$oa1 == m$ea2
m <- m[same | swap]
m[, flip := (ea1 == oa2) & (oa1 == ea2)]
m[, b2h := ifelse(flip, -b2, b2)]
m <- m[!is.na(b1) & !is.na(b2h) & !is.na(se1) & !is.na(se2)]
m[, maf := pmin(eaf1, 1-eaf1)]
m <- m[!is.na(maf) & maf > 0 & maf < 1]
cat("harmonised SNPs for coloc:", nrow(m), "\n")
fwrite(m, file.path(outdir, "coloc_harmonised_input.csv"))

d1 <- list(snp=m$snp, beta=m$b1, varbeta=m$se1^2, type="quant", N=N_EQTL, MAF=m$maf, sdY=1)
d2 <- list(snp=m$snp, beta=m$b2h, varbeta=m$se2^2, type="cc", N=N_OA, s=NCASE/N_OA, MAF=m$maf)
res <- coloc.abf(d1, d2)
cat("\n== COLOC summary ==\n"); print(res$summary)
fwrite(as.data.table(as.list(res$summary)), file.path(outdir, "coloc_VDR_kneeOA_summary.csv"))
pp <- as.data.table(res$results)[order(-SNP.PP.H4)]
fwrite(pp, file.path(outdir, "coloc_VDR_kneeOA_perSNP.csv"))
cat("\ntop SNP by PP.H4:\n"); print(head(pp, 5))

# ---- Steiger directionality ----
cat("\n== Steiger directionality ==\n")
# use the original MR instrument set if available (SuppTable3), else all harmonised cis SNPs
supp <- tryCatch(fread("C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/SuppTable3_VDR_SNPs_Details_and_F_stats.csv"),
                 error=function(e) NULL)
if (!is.null(supp)) {
  setnames(supp, c("SNP","chr","pos","ea","oa","eaf","beta","se","pval","F"))
  ex <- data.frame(SNP=supp$SNP, beta.exposure=as.numeric(supp$beta), se.exposure=as.numeric(supp$se),
                   effect_allele.exposure=toupper(supp$ea), other_allele.exposure=toupper(supp$oa),
                   eaf.exposure=as.numeric(supp$eaf), samplesize.exposure=N_EQTL,
                   units.exposure="SD", id.exposure="VDR_eqtlgen")
  mm <- m[snp %in% supp$SNP]
  ou <- data.frame(SNP=mm$snp, beta.outcome=mm$b2h, se.outcome=mm$se2,
                   effect_allele.outcome=mm$ea1, other_allele.outcome=mm$oa1,
                   eaf.outcome=mm$eaf2, samplesize.outcome=N_OA,
                   ncase.outcome=NCASE, ncontrol.outcome=NCONTROL,
                   units.outcome="log odds", id.outcome="knee_OA")
  dat <- harmonise_data(ex, ou, action=2)
  cat("instruments harmonised:", nrow(dat), " kept:", sum(dat$mr_keep), "\n")
  st <- tryCatch(directionality_test(dat), error=function(e){cat("steiger err:", conditionMessage(e), "\n"); NULL})
  if (!is.null(st)) { print(st); fwrite(as.data.table(st), file.path(outdir, "steiger_VDR_kneeOA.csv")) }
}
cat("\nDONE coloc+steiger\n")
