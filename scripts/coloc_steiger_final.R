suppressMessages({library(data.table); library(coloc); library(TwoSampleMR)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
out <- file.path(base, "outputs")
N_CASE <- 24955; N_CTRL <- 378169; N_OA <- N_CASE + N_CTRL
N_EQTL <- 26609   # OpenGWAS sample size for eqtl-a-ENSG00000111424

eq <- fread(file.path(base, "data/eqtl_region.tsv"))
oa <- fread(file.path(base, "data/oa_region.tsv"))
cat("eQTL rows:", nrow(eq), " OA rows:", nrow(oa), "\n")

eq <- eq[!is.na(rsid) & !is.na(beta) & !is.na(se) & rsid != "."]
oa <- oa[!is.na(rsid) & !is.na(beta) & !is.na(se) & rsid != "."]
eq <- eq[!duplicated(rsid)]; oa <- oa[!duplicated(rsid)]

m <- merge(eq[, .(snp=rsid, ea1=toupper(ea), oa1=toupper(nea), eaf1=as.numeric(eaf),
                  b1=as.numeric(beta), se1=as.numeric(se), p1=as.numeric(p))],
           oa[, .(snp=rsid, ea2=toupper(ea), oa2=toupper(nea), eaf2=as.numeric(eaf),
                  b2=as.numeric(beta), se2=as.numeric(se), p2=as.numeric(p))], by="snp")
m <- m[!(ea1==oa1)]
m <- m[(ea1==ea2 & oa1==oa2) | (ea1==oa2 & oa1==ea2)]
m[, flip := (ea1==oa2) & (oa1==ea2)]
m[, b2h := ifelse(flip, -b2, b2)]
m <- m[!is.na(b1) & !is.na(b2h) & se1>0 & se2>0]
m[, maf := pmin(eaf1, 1-eaf1)]
m <- m[!is.na(maf) & maf>0 & maf<1]
cat("harmonised SNPs for coloc:", nrow(m), "\n")
fwrite(m, file.path(out, "coloc_harmonised_input.csv"))

res <- coloc.abf(list(snp=m$snp, beta=m$b1, varbeta=m$se1^2, type="quant", N=N_EQTL, MAF=m$maf, sdY=1),
                 list(snp=m$snp, beta=m$b2h, varbeta=m$se2^2, type="cc", N=N_OA, s=N_CASE/N_OA, MAF=m$maf))
cat("\n== COLOC summary (VDR eQTL vs knee OA) ==\n"); print(res$summary)
fwrite(as.data.table(as.list(res$summary)), file.path(out, "coloc_VDR_kneeOA_summary.csv"))
pp <- as.data.table(res$results)[order(-SNP.PP.H4)]
fwrite(pp, file.path(out, "coloc_VDR_kneeOA_perSNP.csv"))
cat("\ntop 5 SNPs by PP.H4:\n"); print(head(pp[, .(snp, PP.H4=SNP.PP.H4)][order(-PP.H4)], 5))
cat("\nkey SNP(s):\n"); print(pp[snp %in% c("rs7975232","rs7485057","rs10783222","rs58067800")])

cat("\n== Steiger directionality ==\n")
supp <- tryCatch(fread("C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/SuppTable3_VDR_SNPs_Details_and_F_stats.csv"), error=function(e) NULL)
if (!is.null(supp)) {
  setnames(supp, c("SNP","chr","pos","ea","oa","eaf","beta","se","pval","F"), skip_absent=TRUE)
  ex <- data.frame(SNP=supp$SNP, beta.exposure=as.numeric(supp$beta), se.exposure=as.numeric(supp$se),
                   effect_allele.exposure=toupper(supp$ea), other_allele.exposure=toupper(supp$oa),
                   eaf.exposure=as.numeric(supp$eaf), samplesize.exposure=N_EQTL,
                   units.exposure="SD", id.exposure="VDR_eqtlgen", exposure="VDR expression")
  mm <- m[snp %in% supp$SNP]
  ou <- data.frame(SNP=mm$snp, beta.outcome=mm$b2h, se.outcome=mm$se2,
                   effect_allele.outcome=mm$ea1, other_allele.outcome=mm$oa1,
                   eaf.outcome=mm$eaf2, samplesize.outcome=N_OA, ncase.outcome=N_CASE,
                   ncontrol.outcome=N_CTRL, units.outcome="log odds", id.outcome="knee_OA",
                   outcome="Knee osteoarthritis")
  dat <- harmonise_data(ex, ou, action=2)
  cat("instruments:", nrow(dat), "kept:", sum(dat$mr_keep), "\n")
  st <- tryCatch(directionality_test(dat), error=function(e){cat("steiger err:", conditionMessage(e), "\n"); NULL})
  if (!is.null(st)) { print(st); fwrite(as.data.table(st), file.path(out, "steiger_VDR_kneeOA.csv")) }
}
cat("\nDONE_FINAL\n")
