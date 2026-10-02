# Colocalization + Steiger for VDR -> knee OA using OpenGWAS API (associations endpoint)
suppressMessages({library(ieugwasr); library(TwoSampleMR); library(data.table); library(coloc)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
outdir <- file.path(base, "outputs"); dir.create(outdir, showWarnings=FALSE)
N_CASE <- 24955; N_CTRL <- 378169; N_OA <- N_CASE + N_CTRL

rs <- fread(file.path(base, "data/ens_vdr_region_rsids.txt"), header=FALSE)$V1
cat("region rsIDs:", length(rs), "\n"); flush.console()

cat("== querying VDR cis-eQTL for region SNPs ==\n"); flush.console()
eq <- as.data.table(associations(variants = rs, id = "eqtl-a-ENSG00000111424", proxies = 0))
cat("eQTL rows:", nrow(eq), " cols:", paste(names(eq), collapse=","), "\n"); flush.console()
fwrite(eq, file.path(outdir, "VDR_cis_eQTL_opengwas.csv"))
eq <- eq[!is.na(rsid) & !is.na(beta) & !is.na(se)]
eq <- eq[!duplicated(rsid)]

cat("== querying knee OA for eQTL SNPs ==\n"); flush.console()
oa <- as.data.table(associations(variants = unique(eq$rsid), id = "ebi-a-GCST007090", proxies = 0))
cat("OA rows:", nrow(oa), "\n"); flush.console()
fwrite(oa, file.path(outdir, "kneeOA_region_opengwas.csv"))
oa <- oa[!is.na(rsid) & !is.na(beta) & !is.na(se)]
oa <- oa[!duplicated(rsid)]

e <- data.table(snp=eq$rsid, ea1=toupper(eq$ea), oa1=toupper(eq$nea),
                eaf1=as.numeric(eq$eaf), b1=as.numeric(eq$beta), se1=as.numeric(eq$se),
                n1=as.numeric(eq$n))
o <- data.table(snp=oa$rsid, ea2=toupper(oa$ea), oa2=toupper(oa$nea),
                eaf2=as.numeric(oa$eaf), b2=as.numeric(oa$beta), se2=as.numeric(oa$se))
m <- merge(e, o, by="snp")
m <- m[(ea1==ea2 & oa1==oa2) | (ea1==oa2 & oa1==ea2)]
m <- m[ea1!=oa1]
m[, flip := (ea1==oa2) & (oa1==ea2)]
m[, b2h := ifelse(flip, -b2, b2)]
m <- m[!is.na(b1) & !is.na(b2h) & !is.na(se1) & !is.na(se2) & se1>0 & se2>0]
m[, maf := pmin(eaf1, 1-eaf1)]
m <- m[!is.na(maf) & maf>0 & maf<1]
cat("harmonised SNPs:", nrow(m), "\n"); flush.console()
fwrite(m, file.path(outdir, "coloc_harmonised_input.csv"))

N_EQ <- as.integer(median(m$n1, na.rm=TRUE)); cat("eQTL N:", N_EQ, "\n")
d1 <- list(snp=m$snp, beta=m$b1, varbeta=m$se1^2, type="quant", N=N_EQ, MAF=m$maf, sdY=1)
d2 <- list(snp=m$snp, beta=m$b2h, varbeta=m$se2^2, type="cc", N=N_OA, s=N_CASE/N_OA, MAF=m$maf)
res <- coloc.abf(d1, d2)
cat("\n== COLOC summary ==\n"); print(res$summary)
fwrite(as.data.table(as.list(res$summary)), file.path(outdir, "coloc_VDR_kneeOA_summary.csv"))
pp <- as.data.table(res$results)[order(-SNP.PP.H4)]
fwrite(pp, file.path(outdir, "coloc_VDR_kneeOA_perSNP.csv"))
cat("top SNP by PP.H4:\n"); print(head(pp, 5))

cat("\n== Steiger directionality ==\n"); flush.console()
supp <- tryCatch(fread("C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/SuppTable3_VDR_SNPs_Details_and_F_stats.csv"), error=function(e) NULL)
if (!is.null(supp)) {
  setnames(supp, c("SNP","chr","pos","ea","oa","eaf","beta","se","pval","F"), skip_absent=TRUE)
  ex <- data.frame(SNP=supp$SNP, beta.exposure=as.numeric(supp$beta), se.exposure=as.numeric(supp$se),
                   effect_allele.exposure=toupper(supp$ea), other_allele.exposure=toupper(supp$oa),
                   eaf.exposure=as.numeric(supp$eaf), samplesize.exposure=N_EQ,
                   units.exposure="SD", id.exposure="VDR_eqtlgen")
  mm <- m[snp %in% supp$SNP]
  ou <- data.frame(SNP=mm$snp, beta.outcome=mm$b2h, se.outcome=mm$se2,
                   effect_allele.outcome=mm$ea1, other_allele.outcome=mm$oa1,
                   eaf.outcome=mm$eaf2, samplesize.outcome=N_OA,
                   ncase.outcome=N_CASE, ncontrol.outcome=N_CTRL,
                   units.outcome="log odds", id.outcome="knee_OA")
  dat <- harmonise_data(ex, ou, action=2)
  cat("instruments:", nrow(dat), " kept:", sum(dat$mr_keep), "\n")
  st <- tryCatch(directionality_test(dat), error=function(e){cat("steiger err:", conditionMessage(e), "\n"); NULL})
  if (!is.null(st)) { print(st); fwrite(as.data.table(st), file.path(outdir, "steiger_VDR_kneeOA.csv")) }
}
cat("\nDONE\n")
