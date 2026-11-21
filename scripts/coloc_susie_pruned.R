## coloc.susie on an LD-pruned set (r^2 < 0.8) at the OA peak window
suppressMessages(library(coloc))
suppressMessages(library(susieR))
set.seed(20261121)

DIR <- "C:/Users/liuha/OpenClawVDR/repo/tables/susie_OAwindow"
d <- read.csv(file.path(DIR, "susie_input_pruned.csv"), stringsAsFactors = FALSE)
R <- as.matrix(read.csv(file.path(DIR, "susie_ld_pruned.csv"), row.names = 1, check.names = FALSE))
stopifnot(nrow(R) == nrow(d), identical(rownames(R), d$snp))

cat("variants:", nrow(d), " LD:", paste(dim(R), collapse = "x"), "\n")
cat("max |r| off-diagonal:", round(max(abs(R[lower.tri(R)])), 3), "\n")

N_OA <- 163447 + 1019228
d1 <- list(snp = d$snp, beta = d$b1, varbeta = d$se1^2, N = 31684,
           type = "quant", sdY = 1, LD = R)
d2 <- list(snp = d$snp, beta = d$b2h, varbeta = d$se2^2, N = N_OA,
           type = "cc", s = 163447 / N_OA, LD = R)

abf <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5))$summary
cat("\ncoloc.abf (pruned set):\n"); print(round(abf, 6))

r1 <- runsusie(d1, maxit = 500, verbose = FALSE)
r2 <- runsusie(d2, maxit = 500, verbose = FALSE)

cat("\neQTL converged:", r1$converged, " iterations:", r1$niter, " CS:", length(r1$sets$cs), "\n")
cat("OA   converged:", r2$converged, " iterations:", r2$niter, " CS:", length(r2$sets$cs), "\n")

for (tag in c("eQTL", "OA")) {
  r <- if (tag == "eQTL") r1 else r2
  cat("\n---", tag, "credible sets ---\n")
  if (length(r$sets$cs) == 0) { cat("  none\n"); next }
  for (k in seq_along(r$sets$cs)) {
    idx <- r$sets$cs[[k]]
    lead <- idx[which.max(r$pip[idx])]
    cat(sprintf("  CS%d n=%2d purity=%.3f lead=%-12s pip=%.3f\n",
                k, length(idx), r$sets$purity$min.abs.corr[k], d$snp[lead], r$pip[lead]))
  }
}

if (length(r1$sets$cs) && length(r2$sets$cs)) {
  res <- coloc.susie(r1, r2)
  s <- as.data.frame(res$summary)
  cat("\n=== per-pair PP.H4 (all pairs) ===\n")
  print(s[order(-s$PP.H4.abf), c("idx1", "idx2", "PP.H4.abf")])
  cat(sprintf("\nmax PP.H4 = %.4f ; pairs >0.5: %d ; pairs >0.8: %d\n",
              max(s$PP.H4.abf), sum(s$PP.H4.abf > 0.5), sum(s$PP.H4.abf > 0.8)))
  write.csv(s, file.path(DIR, "coloc_susie_pruned_summary.csv"), row.names = FALSE)
} else {
  cat("\nno credible sets in one trait; susie colocalisation not applicable\n")
}
