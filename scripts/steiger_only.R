suppressMessages({library(ieugwasr); library(TwoSampleMR); library(data.table)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
out <- file.path(base, "outputs")
N_CASE <- 24955; N_CTRL <- 378169; N_OA <- N_CASE + N_CTRL; N_EQTL <- 26609

snps <- c("rs10783222","rs12943633","rs149110519","rs2594836","rs28498283","rs35979828",
          "rs372883","rs4142441","rs58067800","rs7485057","rs7975232")
ex <- as.data.table(associations(variants=snps, id="eqtl-a-ENSG00000111424", proxies=0))
ou <- as.data.table(associations(variants=snps, id="ebi-a-GCST007090", proxies=0))
if ("proxy" %in% names(ex)) ex <- ex[proxy == FALSE | is.na(proxy)]
if ("proxy" %in% names(ou)) ou <- ou[proxy == FALSE | is.na(proxy)]
ex <- ex[!duplicated(rsid)]; ou <- ou[!duplicated(rsid)]
cat("exact exposure SNPs:", nrow(ex), "| exact outcome SNPs:", nrow(ou), "\n")
print(ex[, .(rsid, ea, nea, eaf, beta, se, p)])
fwrite(ex, file.path(out, "steiger_exposure_snps.csv"))
fwrite(ou, file.path(out, "steiger_outcome_snps.csv"))

e <- data.frame(SNP=ex$rsid, beta.exposure=as.numeric(ex$beta), se.exposure=as.numeric(ex$se),
                effect_allele.exposure=toupper(ex$ea), other_allele.exposure=toupper(ex$nea),
                eaf.exposure=as.numeric(ex$eaf), pval.exposure=as.numeric(ex$p),
                samplesize.exposure=N_EQTL, units.exposure="SD",
                id.exposure="VDR_eqtlgen", exposure="VDR expression")
o <- data.frame(SNP=ou$rsid, beta.outcome=as.numeric(ou$beta), se.outcome=as.numeric(ou$se),
                effect_allele.outcome=toupper(ou$ea), other_allele.outcome=toupper(ou$nea),
                eaf.outcome=as.numeric(ou$eaf), pval.outcome=as.numeric(ou$p),
                samplesize.outcome=N_OA, ncase.outcome=N_CASE, ncontrol.outcome=N_CTRL,
                units.outcome="log odds", id.outcome="knee_OA", outcome="Knee osteoarthritis")
dat <- suppressMessages(harmonise_data(e, o, action=2))
cat("harmonised instruments:", nrow(dat), "| mr_keep:", sum(dat$mr_keep), "\n")
print(dat[, c("SNP","beta.exposure","beta.outcome","mr_keep")])
st <- tryCatch(directionality_test(dat), error=function(err){cat("steiger err:", conditionMessage(err), "\n"); NULL})
if (!is.null(st)) { print(st); fwrite(as.data.table(st), file.path(out, "steiger_VDR_kneeOA.csv")) }

# binary-correct r.outcome (case-control) and re-run
library(TwoSampleMR)
dat$r.exposure <- get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome, dat$ncase.outcome, dat$ncontrol.outcome, 0.15)
st2 <- tryCatch(directionality_test(dat), error=function(err){cat("steiger2 err:", conditionMessage(err), "\n"); NULL})
cat("\n-- Steiger with binary-corrected r.outcome --\n")
if (!is.null(st2)) { print(st2); fwrite(as.data.table(st2), file.path(out, "steiger_VDR_kneeOA_binarycorrected.csv")) }
cat("\nDONE_STEIGER\n")
