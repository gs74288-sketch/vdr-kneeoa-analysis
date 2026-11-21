## Inspect the coloc.susie result: which variants are in the shared credible sets?
suppressMessages(library(coloc))
suppressMessages(library(susieR))

DIR <- "C:/Users/liuha/OpenClawVDR/repo/tables/susie_OAwindow"
d <- read.csv(file.path(DIR, "susie_input.csv"), stringsAsFactors = FALSE)
fit <- readRDS(file.path(DIR, "coloc_susie_full.rds"))
r1 <- fit$r1; r2 <- fit$r2; res <- fit$res

cat("=== credible set counts ===\n")
cat("eQTL:", length(r1$sets$cs), " OA:", length(r2$sets$cs), "\n")

show <- function(r, tag) {
  cat(sprintf("\n--- %s credible sets ---\n", tag))
  for (k in seq_along(r$sets$cs)) {
    idx <- r$sets$cs[[k]]
    lead <- idx[which.max(r$pip[idx])]
    cat(sprintf("CS%-2d  n=%3d  purity=%.3f  lead=%s  leadPIP=%.3f  maxPIP=%.3f\n",
                k, length(idx), r$sets$purity$min.abs.corr[k],
                d$snp[lead], r$pip[lead], max(r$pip[idx])))
  }
}
show(r1, "eQTL")
show(r2, "OA")

cat("\n=== top 10 PIP, eQTL ===\n")
o <- order(-r1$pip)[1:10]
print(data.frame(snp = d$snp[o], pip = round(r1$pip[o], 3),
                 p_eqtl = signif(d$p1[o], 2), p_oa = signif(d$p2[o], 2)))

cat("\n=== top 10 PIP, OA ===\n")
o2 <- order(-r2$pip)[1:10]
print(data.frame(snp = d$snp[o2], pip = round(r2$pip[o2], 3),
                 p_oa = signif(d$p2[o2], 2), p_eqtl = signif(d$p1[o2], 2)))

s <- as.data.frame(res$summary)
cat("\n=== pairs with PP.H4 > 0.5 ===\n")
top <- s[s$PP.H4.abf > 0.5, ]
for (i in seq_len(nrow(top))) {
  i1 <- top$idx1[i]; i2 <- top$idx2[i]
  cs1 <- r1$sets$cs[[i1]]; cs2 <- r2$sets$cs[[i2]]
  cat(sprintf("\npair PP.H4=%.4f  eQTL CS%d (n=%d)  OA CS%d (n=%d)\n",
              top$PP.H4.abf[i], i1, length(cs1), i2, length(cs2)))
  cat("  eQTL CS leads:", paste(d$snp[cs1[order(-r1$pip[cs1])][1:min(3,length(cs1))]], collapse=", "), "\n")
  cat("  OA   CS leads:", paste(d$snp[cs2[order(-r2$pip[cs2])][1:min(3,length(cs2))]], collapse=", "), "\n")
  cat("  shared variants:", paste(d$snp[intersect(cs1, cs2)], collapse=", "), "\n")
  cat("  eQTL p of OA lead:", signif(d$p1[cs2[which.max(r2$pip[cs2])]], 2), "\n")
  cat("  OA   p of eQTL lead:", signif(d$p2[cs1[which.max(r1$pip[cs1])]], 3), "\n")
}

cat("\n=== sanity: pi and niter ===\n")
cat("eQTL niter:", r1$niter, " L:", r1$L, "  OA niter:", r2$niter, " L:", r2$L, "\n")
cat("eQTL sigma2:", round(r1$sigma2, 2), " OA sigma2:", round(r2$sigma2, 2), "\n")
