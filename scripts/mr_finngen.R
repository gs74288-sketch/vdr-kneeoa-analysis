suppressMessages({library(data.table); library(TwoSampleMR)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
fin  <- "C:/Users/liuha/OpenClawVDR/finngen"
out  <- file.path(base, "outputs")

run_one <- function(f, label) {
  if (!file.exists(f)) { cat("missing:", f, "\n"); return(NULL) }
  ff <- fread(f)
  setnames(ff, gsub("^#", "", names(ff)))
  cat("\n=====", label, "| rows:", nrow(ff), "=====\n")
  # exact rsID matching against the instrument list
  supp <- fread("C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/SuppTable3_VDR_SNPs_Details_and_F_stats.csv")
  setnames(supp, c("SNP","chr","pos","ea","oa","eaf","beta","se","pval","F"), skip_absent=TRUE)
  want <- supp$SNP
  ff[, snp_exact := sapply(rsids, function(x) { k <- intersect(strsplit(x, ",")[[1]], want); if (length(k)) k[1] else NA_character_ })]
  ff <- ff[!is.na(snp_exact)]
  ff <- ff[!duplicated(snp_exact)]
  cat("matched instruments:", nrow(ff), ":", paste(ff$snp_exact, collapse=","), "\n")
  if (nrow(ff) == 0) return(NULL)

  ou <- data.frame(SNP=ff$snp_exact, beta.outcome=as.numeric(ff$beta), se.outcome=as.numeric(ff$sebeta),
                   effect_allele.outcome=toupper(ff$alt), other_allele.outcome=toupper(ff$ref),
                   eaf.outcome=as.numeric(ff$af_alt), pval.outcome=as.numeric(ff$pval),
                   samplesize.outcome=NA, units.outcome="log odds",
                   id.outcome=label, outcome=label)
  ex <- data.frame(SNP=supp$SNP, beta.exposure=as.numeric(supp$beta), se.exposure=as.numeric(supp$se),
                   effect_allele.exposure=toupper(supp$ea), other_allele.exposure=toupper(supp$oa),
                   eaf.exposure=as.numeric(supp$eaf), pval.exposure=as.numeric(supp$pval),
                   samplesize.exposure=31684, units.exposure="SD",
                   id.exposure="VDR_eqtlgen", exposure="VDR expression")
  dat <- as.data.frame(suppressMessages(harmonise_data(ex, ou, action=2)))
  dat <- dat[dat$mr_keep, , drop=FALSE]
  cat("harmonised:", nrow(dat), "\n")
  if (nrow(dat) < 1) return(NULL)
  fwrite(dat, file.path(out, paste0("finngen_", label, "_harmonised.csv")))
  ml <- if (nrow(dat) >= 2) c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode") else "mr_wald_ratio"
  res <- generate_odds_ratios(mr(dat, method_list=ml))
  print(res[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  fwrite(res, file.path(out, paste0("finngen_", label, "_results.csv")))
  if (nrow(dat) >= 2) { print(mr_heterogeneity(dat)); print(mr_pleiotropy_test(dat)) }
  invisible(res)
}

r1 <- run_one(file.path(fin, "M13_ARTHROSIS_KNEE_instrument_snps.tsv"), "FinnGen_knee_arthrosis")
r2 <- run_one(file.path(fin, "M13_ARTHROSIS_KNEE_PRIM_KNEESURG_instrument_snps.tsv"), "FinnGen_knee_arthrosis_kneesurgery")
cat("\nDONE_FINNGEN_MR\n")
