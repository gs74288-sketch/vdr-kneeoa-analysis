suppressMessages({library(data.table); library(TwoSampleMR)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
out  <- file.path(base, "outputs")
fin  <- "C:/Users/liuha/OpenClawVDR/finngen"

ins <- fread(file.path(out, "cis_ldclumped_instruments.csv"))
ex <- data.frame(SNP=ins$SNP, beta.exposure=as.numeric(ins$beta), se.exposure=as.numeric(ins$se),
                 effect_allele.exposure=toupper(ins$ea), other_allele.exposure=toupper(ins$nea),
                 eaf.exposure=as.numeric(ins$eaf), pval.exposure=as.numeric(ins$p),
                 samplesize.exposure=31684, units.exposure="SD",
                 id.exposure="VDR_cis_LDclumped", exposure="VDR expression (cis, LD-clumped)")

run <- function(f, lab) {
  if (!file.exists(f)) { cat("missing", f, "\n"); return(NULL) }
  ff <- fread(f); setnames(ff, gsub("^#", "", names(ff)))
  ff <- ff[rsids %in% ins$SNP]
  ff <- ff[!duplicated(rsids)]
  cat("\n=====", lab, "| matched:", nrow(ff), "=====\n")
  if (!nrow(ff)) return(NULL)
  ou <- data.frame(SNP=ff$rsids, beta.outcome=as.numeric(ff$beta), se.outcome=as.numeric(ff$sebeta),
                   effect_allele.outcome=toupper(ff$alt), other_allele.outcome=toupper(ff$ref),
                   eaf.outcome=as.numeric(ff$af_alt), pval.outcome=as.numeric(ff$pval),
                   units.outcome="log odds", id.outcome=lab, outcome=lab)
  dat <- as.data.frame(suppressMessages(harmonise_data(ex, ou, action=2)))
  dat <- dat[dat$mr_keep, , drop=FALSE]
  cat("harmonised:", nrow(dat), "\n")
  if (!nrow(dat)) return(NULL)
  fwrite(dat, file.path(out, paste0("cisMR_LDclumped_", lab, "_harmonised.csv")))
  ml <- if (nrow(dat) >= 2) c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode") else "mr_wald_ratio"
  res <- generate_odds_ratios(mr(dat, method_list=ml))
  print(res[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  fwrite(res, file.path(out, paste0("cisMR_LDclumped_", lab, "_results.csv")))
  if (nrow(dat) >= 2) print(mr_heterogeneity(dat))
  invisible(res)
}

run(file.path(fin, "M13_ARTHROSIS_KNEE_cisld_snps.tsv"), "FinnGen_knee_arthrosis")
run(file.path(fin, "M13_ARTHROSIS_KNEE_PRIM_KNEESURG_cisld_snps.tsv"), "FinnGen_knee_arthrosis_kneesurgery")
cat("\nDONE_CIS_LDCLUMPED_FINNGEN\n")
