## Smoke test: can we pull a 1000G EUR LD matrix through OpenGWAS?
suppressMessages(library(ieugwasr))

snps <- c("rs7975232", "rs11168351", "rs4506714", "rs61919782", "rs7307566")
cat("testing", length(snps), "SNPs\n")

m <- try(ld_matrix(snps, pop = "EUR", with_alleles = FALSE), silent = TRUE)
if (inherits(m, "try-error")) {
  cat("FAILED\n")
  cat(as.character(m))
} else {
  cat("OK  dim =", paste(dim(m), collapse = " x "), "\n")
  cat("rownames head:\n"); print(head(rownames(m), 3))
  cat("diag:\n"); print(round(diag(m), 3))
  cat("off-diagonal range:", round(range(m[lower.tri(m)]), 3), "\n")
}
