## coloc.susie at the OA peak window (chr12:47.91-48.11 Mb, GRCh38)
## LD: 1000G phase3 EUR, 503 samples, estimated from regional genotypes.
suppressMessages(library(coloc))
suppressMessages(library(susieR))

set.seed(20261121)

DIR <- "C:/Users/liuha/OpenClawVDR/repo/tables/susie_OAwindow"
d <- read.csv(file.path(DIR, "susie_input.csv"), stringsAsFactors = FALSE)
R <- as.matrix(read.csv(file.path(DIR, "susie_ld.csv"), row.names = 1, check.names = FALSE))
snps <- readLines(file.path(DIR, "susie_ld_snps.txt"))

cat("variants:", nrow(d), " LD:", paste(dim(R), collapse = "x"), "\n")
stopifnot(identical(d$snp, snps), nrow(R) == nrow(d))

N_EQTL <- 31684
N_CASE <- 163447
N_CTRL <- 1019228
N_OA <- N_CASE + N_CTRL
S_OA <- N_CASE / N_OA

d1 <- list(snp = d$snp, beta = d$b1, varbeta = d$se1^2, N = N_EQTL,
           type = "quant", sdY = 1, LD = R)
d2 <- list(snp = d$snp, beta = d$b2h, varbeta = d$se2^2, N = N_OA,
           type = "cc", s = S_OA, LD = R)

## reference point: single-causal-variant model on exactly the same variants
abf <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5))$summary
cat("\n=== coloc.abf on the matched 701 variants ===\n")
print(round(abf, 6))

## susie fine-mapping per trait
r1 <- runsusie(d1, maxit = 1000, verbose = FALSE)
r2 <- runsusie(d2, maxit = 1000, verbose = FALSE)

summ <- function(r, tag) {
  cat(sprintf("\n%-6s credible sets: %d\n", tag, length(r$sets$cs)))
  if (length(r$sets$cs)) {
    for (k in seq_along(r$sets$cs)) {
      idx <- r$sets$cs[[k]]
      lead <- idx[which.max(r$pip[idx])]
      cat(sprintf("   CS%d  n=%3d  purity=%.3f  lead=%s  pip=%.3f\n",
                  k, length(idx), r$sets$purity$min.abs.corr[k],
                  d$snp[lead], r$pip[lead]))
    }
  }
}
summ(r1, "eQTL")
summ(r2, "OA")

res <- try(coloc.susie(r1, r2), silent = TRUE)
if (inherits(res, "try-error")) {
  cat("\ncoloc.susie FAILED:\n"); cat(as.character(res)); quit(status = 1)
}

s <- res$summary
cat("\n=== coloc.susie: per-pair results ===\n")
print(s)

if (length(s) && !is.null(s$PP.H4.abf)) {
  best <- which.max(s$PP.H4.abf)
  cat(sprintf("\nmax PP.H4 across pairs = %.4f", max(s$PP.H4.abf)))
  cat(sprintf("  (pair %d: idx1=%s idx2=%s)\n", best, s$idx1[best], s$idx2[best]))
  cat(sprintf("pairs with PP.H4 > 0.5 : %d\n", sum(s$PP.H4.abf > 0.5)))
  cat(sprintf("pairs with PP.H4 > 0.8 : %d\n", sum(s$PP.H4.abf > 0.8)))
}

write.csv(as.data.frame(s), file.path(DIR, "coloc_susie_summary.csv"), row.names = FALSE)
saveRDS(list(r1 = r1, r2 = r2, res = res), file.path(DIR, "coloc_susie_full.rds"))
cat("\nwritten coloc_susie_summary.csv\n")
