suppressMessages({library(ieugwasr); library(TwoSampleMR); library(data.table); library(coloc)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
out <- file.path(base, "outputs")
N_CASE <- 24955; N_CTRL <- 378169; N_OA <- N_CASE + N_CTRL
N_EQTL <- 31684   # eQTLGen whole-blood cis-eQTL sample size (Vosa et al. 2021)

eq <- fread(file.path(base, "data/eqtl_region.tsv"))
oa <- fread(file.path(base, "data/oa_region.tsv"))
eq <- eq[!is.na(rsid) & rsid != "." & !is.na(beta) & !is.na(se) & se > 0]
oa <- oa[!is.na(rsid) & rsid != "." & !is.na(beta) & !is.na(se) & se > 0]
eq <- eq[!duplicated(rsid)]; oa <- oa[!duplicated(rsid)]

cat("=== 1) PURE-CIS instrument MR (chr12:46.8-48.9Mb) ===\n")
cand <- eq[p < 5e-5]
cand[, F := (beta/se)^2]
cand <- cand[F > 10]
cat("candidate cis SNPs (P<5e-5, F>10):", nrow(cand), "\n")
cl <- tryCatch(ld_clump(data.frame(rsid=cand$rsid, p=cand$p, id="eqtl-a-ENSG00000111424"),
                        clump_kb=10000, clump_r2=0.001, pop="EUR"),
               error=function(e){cat("ld_clump err:", conditionMessage(e), "\n"); NULL})
if (!is.null(cl)) {
  cat("cis instruments after clumping:", nrow(cl), ":", paste(cl$rsid, collapse=","), "\n")
  ex <- data.frame(SNP=cl$rsid, beta.exposure=cand[match(cl$rsid, rsid)]$beta,
                   se.exposure=cand[match(cl$rsid, rsid)]$se,
                   effect_allele.exposure=toupper(cand[match(cl$rsid, rsid)]$ea),
                   other_allele.exposure=toupper(cand[match(cl$rsid, rsid)]$nea),
                   eaf.exposure=cand[match(cl$rsid, rsid)]$eaf, pval.exposure=cand[match(cl$rsid, rsid)]$p,
                   samplesize.exposure=N_EQTL, units.exposure="SD",
                   id.exposure="VDR_cis", exposure="VDR expression (cis)")
  ou <- oa[rsid %in% cl$rsid]
  ou <- data.frame(SNP=ou$rsid, beta.outcome=ou$beta, se.outcome=ou$se,
                   effect_allele.outcome=toupper(ou$ea), other_allele.outcome=toupper(ou$nea),
                   eaf.outcome=ou$eaf, pval.outcome=ou$p, samplesize.outcome=N_OA,
                   ncase.outcome=N_CASE, ncontrol.outcome=N_CTRL, units.outcome="log odds",
                   id.outcome="knee_OA", outcome="Knee osteoarthritis")
  dat <- suppressMessages(harmonise_data(ex, ou, action=2))
  dat <- dat[mr_keep]
  cat("harmonised cis instruments:", nrow(dat), "\n")
  fwrite(dat, file.path(out, "cisMR_harmonised.csv"))
  res <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  or <- generate_odds_ratios(res)
  print(or[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  fwrite(or, file.path(out, "cisMR_results.csv"))
  het <- mr_heterogeneity(dat); print(het)
  pl <- mr_pleiotropy_test(dat); print(pl)
  fwrite(as.data.table(het), file.path(out, "cisMR_heterogeneity.csv"))
  fwrite(as.data.table(pl), file.path(out, "cisMR_pleiotropy.csv"))
  loo <- mr_leaveoneout(dat); fwrite(as.data.table(loo), file.path(out, "cisMR_leaveoneout.csv"))
}

cat("\n=== 2) coloc with N_EQTL = 31684 ===\n")
m <- merge(eq[, .(snp=rsid, ea1=toupper(ea), oa1=toupper(nea), eaf1=as.numeric(eaf),
                  b1=as.numeric(beta), se1=as.numeric(se))],
           oa[, .(snp=rsid, ea2=toupper(ea), oa2=toupper(nea), eaf2=as.numeric(eaf),
                  b2=as.numeric(beta), se2=as.numeric(se))], by="snp")
m <- m[ea1 != oa1]
m <- m[(ea1==ea2 & oa1==oa2) | (ea1==oa2 & oa1==ea2)]
m[, flip := (ea1==oa2) & (oa1==ea2)]
m[, b2h := ifelse(flip, -b2, b2)]
m[, maf := pmin(eaf1, 1-eaf1)]
m <- m[!is.na(maf) & maf>0 & maf<1 & se1>0 & se2>0]
cat("coloc SNPs:", nrow(m), "\n")
res <- coloc.abf(list(snp=m$snp, beta=m$b1, varbeta=m$se1^2, type="quant", N=N_EQTL, MAF=m$maf, sdY=1),
                 list(snp=m$snp, beta=m$b2h, varbeta=m$se2^2, type="cc", N=N_OA, s=N_CASE/N_OA, MAF=m$maf))
print(res$summary)
fwrite(as.data.table(as.list(res$summary)), file.path(out, "coloc_VDR_kneeOA_summary_N31684.csv"))

cat("\n=== 3) Steiger with N_EQTL = 31684 ===\n")
supp <- fread("C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/SuppTable3_VDR_SNPs_Details_and_F_stats.csv")
setnames(supp, c("SNP","chr","pos","ea","oa","eaf","beta","se","pval","F"), skip_absent=TRUE)
snps <- supp$SNP
ex2 <- as.data.table(associations(variants=snps, id="eqtl-a-ENSG00000111424", proxies=0))
ou2 <- as.data.table(associations(variants=snps, id="ebi-a-GCST007090", proxies=0))
if ("proxy" %in% names(ex2)) ex2 <- ex2[proxy==FALSE | is.na(proxy)]
if ("proxy" %in% names(ou2)) ou2 <- ou2[proxy==FALSE | is.na(proxy)]
E <- data.frame(SNP=ex2$rsid, beta.exposure=ex2$beta, se.exposure=ex2$se,
                effect_allele.exposure=toupper(ex2$ea), other_allele.exposure=toupper(ex2$nea),
                eaf.exposure=ex2$eaf, pval.exposure=ex2$p, samplesize.exposure=N_EQTL,
                units.exposure="SD", id.exposure="VDR_eqtlgen", exposure="VDR expression")
O <- data.frame(SNP=ou2$rsid, beta.outcome=ou2$beta, se.outcome=ou2$se,
                effect_allele.outcome=toupper(ou2$ea), other_allele.outcome=toupper(ou2$nea),
                eaf.outcome=ou2$eaf, pval.outcome=ou2$p, samplesize.outcome=N_OA,
                ncase.outcome=N_CASE, ncontrol.outcome=N_CTRL, units.outcome="log odds",
                id.outcome="knee_OA", outcome="Knee osteoarthritis")
dat2 <- suppressMessages(harmonise_data(E, O, action=2))
dat2$r.exposure <- get_r_from_bsen(dat2$beta.exposure, dat2$se.exposure, dat2$samplesize.exposure)
dat2$r.outcome <- get_r_from_lor(dat2$beta.outcome, dat2$eaf.outcome, dat2$ncase.outcome, dat2$ncontrol.outcome, 0.15)
st <- directionality_test(dat2); print(st)
fwrite(as.data.table(st), file.path(out, "steiger_VDR_kneeOA_N31684.csv"))
cat("\nDONE_CIS_RECOLOC\n")
