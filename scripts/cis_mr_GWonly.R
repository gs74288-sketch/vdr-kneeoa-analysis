# Sensitivity analysis: cis-restricted IVW restricted to instruments that reach
# genome-wide significance for VDR expression (P < 5e-8).
#
# Input : tables/cisMR_LDclumped_discovery_harmonised.csv  (4 LD-clumped cis variants)
# Output: tables/cisMR_GWonly_harmonised.csv
#         tables/cisMR_GWonly_results.csv
#
# Both fixed-effect (mr_ivw_fe) and multiplicative random-effects (mr_ivw_mre)
# estimates are written, because with three variants the heterogeneity model
# changes the width of the confidence interval but not the point estimate.
#
# TwoSampleMR 0.7.9: mr_ivw() has no `model` argument; use mr_ivw_fe / mr_ivw_mre.

suppressMessages({library(data.table); library(TwoSampleMR)})

harm <- fread("tables/cisMR_LDclumped_discovery_harmonised.csv")
harm <- harm[mr_keep == TRUE]
cat("LD-clumped cis variants:", nrow(harm), "\n")

gw <- harm[pval.exposure < 5e-8]
cat("reaching genome-wide significance (P < 5e-8):", nrow(gw), "\n")
cat("excluded:", paste(setdiff(harm$SNP, gw$SNP), collapse = ", "), "\n")

fwrite(gw, "tables/cisMR_GWonly_harmonised.csv")

fe  <- mr_ivw_fe(gw$beta.exposure,  gw$beta.outcome,  gw$se.exposure,  gw$se.outcome)
mre <- mr_ivw_mre(gw$beta.exposure, gw$beta.outcome, gw$se.exposure, gw$se.outcome)

out <- rbind(
  data.table(
    id.exposure = "VDR_cis_GWonly", id.outcome = "GCST007090",
    outcome = "Knee OA (discovery)",
    exposure = "VDR expression (cis, LD-clumped, genome-wide significant)",
    method = "Inverse variance weighted (random effects)", nsnp = nrow(gw),
    b = mre$b, se = mre$se, pval = mre$pval,
    lo_ci = mre$b - 1.96 * mre$se, up_ci = mre$b + 1.96 * mre$se,
    or = exp(mre$b), or_lci95 = exp(mre$b - 1.96 * mre$se),
    or_uci95 = exp(mre$b + 1.96 * mre$se)),
  data.table(
    id.exposure = "VDR_cis_GWonly", id.outcome = "GCST007090",
    outcome = "Knee OA (discovery)",
    exposure = "VDR expression (cis, LD-clumped, genome-wide significant)",
    method = "Inverse variance weighted (fixed effect)", nsnp = nrow(gw),
    b = fe$b, se = fe$se, pval = fe$pval,
    lo_ci = fe$b - 1.96 * fe$se, up_ci = fe$b + 1.96 * fe$se,
    or = exp(fe$b), or_lci95 = exp(fe$b - 1.96 * fe$se),
    or_uci95 = exp(fe$b + 1.96 * fe$se))
)

fwrite(out, "tables/cisMR_GWonly_results.csv")

cat(sprintf("\nrandom effects: OR %.3f (95%% CI %.3f-%.3f), P = %.3f\n",
            exp(mre$b), exp(mre$b - 1.96*mre$se), exp(mre$b + 1.96*mre$se), mre$pval))
cat(sprintf("fixed effect  : OR %.3f (95%% CI %.3f-%.3f), P = %.3f\n",
            exp(fe$b), exp(fe$b - 1.96*fe$se), exp(fe$b + 1.96*fe$se), fe$pval))
cat(sprintf("Cochran Q = %.2f on %d df, P = %.3f\n", fe$Q, fe$Q_df, fe$Q_pval))
