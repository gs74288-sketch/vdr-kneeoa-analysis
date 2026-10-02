suppressMessages({library(data.table); library(TwoSampleMR)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
out <- file.path(base, "outputs")
N_CASE <- 24955; N_CTRL <- 378169; N_OA <- N_CASE + N_CTRL; N_EQTL <- 31684

eq <- fread(file.path(base, "data/eqtl_region.tsv"))
oa <- fread(file.path(base, "data/oa_region.tsv"))
eq <- eq[!is.na(rsid) & rsid != "." & !is.na(beta) & !is.na(se) & se > 0 & !is.na(p)]
oa <- oa[!is.na(rsid) & rsid != "." & !is.na(beta) & !is.na(se) & se > 0]
eq <- eq[!duplicated(rsid)]; oa <- oa[!duplicated(rsid)]
cand <- eq[p < 5e-5]; cand[, F := (beta/se)^2]; cand <- cand[F > 10 & !is.na(position)]
setorder(cand, p)
cat("candidate cis SNPs:", nrow(cand), "\n")

clump_local <- function(dt, kb){
  selp <- c(); keep <- logical(nrow(dt))
  for (i in seq_len(nrow(dt))) {
    if (all(abs(dt$position[i] - selp) > kb*1000)) { keep[i] <- TRUE; selp <- c(selp, dt$position[i]) }
  }
  dt[keep]
}

run_mr <- function(ins, label){
  ex <- data.frame(SNP=ins$rsid, beta.exposure=ins$beta, se.exposure=ins$se,
                   effect_allele.exposure=toupper(ins$ea), other_allele.exposure=toupper(ins$nea),
                   eaf.exposure=ins$eaf, pval.exposure=ins$p, samplesize.exposure=N_EQTL,
                   units.exposure="SD", id.exposure="VDR_cis", exposure="VDR expression (cis)")
  oo <- oa[rsid %in% ins$rsid]
  ou <- data.frame(SNP=oo$rsid, beta.outcome=oo$beta, se.outcome=oo$se,
                   effect_allele.outcome=toupper(oo$ea), other_allele.outcome=toupper(oo$nea),
                   eaf.outcome=oo$eaf, pval.outcome=oo$p, samplesize.outcome=N_OA,
                   ncase.outcome=N_CASE, ncontrol.outcome=N_CTRL, units.outcome="log odds",
                   id.outcome="knee_OA", outcome="Knee osteoarthritis")
  dat <- as.data.frame(suppressMessages(harmonise_data(ex, ou, action=2)))
  dat <- dat[dat$mr_keep, , drop=FALSE]
  if (nrow(dat) == 0) { cat(label, ": no harmonised instruments\n"); return(NULL) }
  cat("\n--", label, "-- instruments:", nrow(dat), ":", paste(dat$SNP, collapse=","), "\n")
  ml <- if (nrow(dat) >= 2) c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode") else "mr_wald_ratio"
  res <- generate_odds_ratios(mr(dat, method_list=ml))
  print(res[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  fwrite(res, file.path(out, paste0("cisMR_", label, "_results.csv")))
  fwrite(dat, file.path(out, paste0("cisMR_", label, "_harmonised.csv")))
  if (nrow(dat) >= 2) {
    print(mr_heterogeneity(dat)); print(mr_pleiotropy_test(dat))
  }
  invisible(res)
}

ins10 <- clump_local(cand, 10)
ins100 <- clump_local(cand, 100)
cat("distance-clumped instruments: 10kb =", nrow(ins10), "; 100kb =", nrow(ins100), "\n")
r10 <- run_mr(ins10, "10kb")
r100 <- run_mr(ins100, "100kb")
top <- cand[1]
rtop <- run_mr(top, "topSNP")
cat("\nDONE_CIS_LOCAL\n")
