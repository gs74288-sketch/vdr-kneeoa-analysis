suppressMessages({library(data.table); library(TwoSampleMR); library(ieugwasr)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
out  <- file.path(base, "outputs")

ins <- fread(file.path(out, "cis_ldclumped_instruments.csv"))
cat("LD-clumped cis instruments:", nrow(ins), "\n"); print(ins)

ex <- data.frame(SNP=ins$SNP, beta.exposure=as.numeric(ins$beta), se.exposure=as.numeric(ins$se),
                 effect_allele.exposure=toupper(ins$ea), other_allele.exposure=toupper(ins$nea),
                 eaf.exposure=as.numeric(ins$eaf), pval.exposure=as.numeric(ins$p),
                 samplesize.exposure=31684, units.exposure="SD",
                 id.exposure="VDR_cis_LDclumped", exposure="VDR expression (cis, LD-clumped)")

# ---- discovery outcome ----
oa <- as.data.table(associations(variants=ins$SNP, id="ebi-a-GCST007090", proxies=0))
if ("proxy" %in% names(oa)) oa <- oa[proxy == FALSE | is.na(proxy)]
oa <- oa[!duplicated(rsid)]
cat("\ndiscovery outcome rows:", nrow(oa), "\n"); print(oa[, .(rsid, ea, nea, eaf, beta, se, p)])
ou <- data.frame(SNP=oa$rsid, beta.outcome=as.numeric(oa$beta), se.outcome=as.numeric(oa$se),
                 effect_allele.outcome=toupper(oa$ea), other_allele.outcome=toupper(oa$nea),
                 eaf.outcome=as.numeric(oa$eaf), pval.outcome=as.numeric(oa$p),
                 samplesize.outcome=403124, ncase.outcome=24955, ncontrol.outcome=378169,
                 units.outcome="log odds", id.outcome="GCST007090", outcome="Knee OA (discovery)")

dat <- as.data.frame(suppressMessages(harmonise_data(ex, ou, action=2)))
dat <- dat[dat$mr_keep, , drop=FALSE]
cat("\nharmonised:", nrow(dat), "\n"); print(dat[, c("SNP","beta.exposure","beta.outcome")])
fwrite(dat, file.path(out, "cisMR_LDclumped_discovery_harmonised.csv"))
ml <- if (nrow(dat) >= 2) c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode") else "mr_wald_ratio"
res <- generate_odds_ratios(mr(dat, method_list=ml))
print(res[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
fwrite(res, file.path(out, "cisMR_LDclumped_discovery_results.csv"))
if (nrow(dat) >= 2) { print(mr_heterogeneity(dat)); print(mr_pleiotropy_test(dat)) }
cat("\nDONE_CIS_LDCLUMPED\n")
