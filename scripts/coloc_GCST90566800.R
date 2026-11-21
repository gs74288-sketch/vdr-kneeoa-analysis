## Colocalisation of VDR cis-eQTL (eQTLGen) with knee OA (GCST90566800)
## dataset1 = VDR expression (quant, eQTLGen n = 31,684)
## dataset2 = knee OA (cc, 2025 multi-ancestry GWAS, European-dominant)

suppressMessages(library(coloc))

IN  <- "C:/Users/liuha/OpenClawVDR/repo/tables/coloc_harmonised_input_GCST90566800.csv"
OUT <- "C:/Users/liuha/OpenClawVDR/repo/tables"

d <- read.csv(IN, stringsAsFactors = FALSE)
cat("SNPs:", nrow(d), "\n")

N_EQTL <- 31684
N_CASE <- median(d$ncase, na.rm = TRUE)
N_CTRL <- median(d$ncontrol, na.rm = TRUE)
N_OA   <- N_CASE + N_CTRL
S_OA   <- N_CASE / N_OA
cat(sprintf("OA N = %.0f (cases %.0f, controls %.0f), s = %.4f\n", N_OA, N_CASE, N_CTRL, S_OA))

d1 <- list(snp = d$snp, beta = d$b1, varbeta = d$se1^2, N = N_EQTL,
           type = "quant", sdY = 1, MAF = d$maf)
d2 <- list(snp = d$snp, beta = d$b2h, varbeta = d$se2^2, N = N_OA,
           type = "cc", s = S_OA, MAF = d$maf)

run <- function(p12, tag) {
  r <- suppressMessages(coloc.abf(d1, d2, p12 = p12))
  s <- r$summary
  cat(sprintf("\n=== %s (p12 = %g) ===\n", tag, p12))
  cat(sprintf("  nsnps  %d\n", s[["nsnps"]]))
  cat(sprintf("  PP.H0  %.4g\n  PP.H1  %.4g\n  PP.H2  %.4g\n  PP.H3  %.4g\n  PP.H4  %.4g\n",
              s[["PP.H0.abf"]], s[["PP.H1.abf"]], s[["PP.H2.abf"]],
              s[["PP.H3.abf"]], s[["PP.H4.abf"]]))
  res <- r$results
  top <- res[which.max(res$SNP.PP.H4), ]
  cat(sprintf("  top shared candidate: %s  SNP.PP.H4 = %.3f\n", top$snp, top$SNP.PP.H4))
  list(summary = s, results = res)
}

a <- run(1e-5, "primary")
b <- run(1e-6, "sensitivity p12=1e-6")
c <- run(1e-4, "sensitivity p12=1e-4")

write.csv(as.data.frame(t(a$summary)), file.path(OUT, "coloc_GCST90566800_summary.csv"),
          row.names = TRUE)
write.csv(a$results, file.path(OUT, "coloc_GCST90566800_perSNP.csv"), row.names = FALSE)

sens <- rbind(
  data.frame(p12 = 1e-6, t(as.data.frame(b$summary))),
  data.frame(p12 = 1e-5, t(as.data.frame(a$summary))),
  data.frame(p12 = 1e-4, t(as.data.frame(c$summary)))
)
write.csv(sens, file.path(OUT, "coloc_GCST90566800_p12_sensitivity.csv"), row.names = FALSE)

## Where do the two signals peak?
cat("\n=== signal peaks ===\n")
cat(sprintf("eQTL  min p: %s (%.3g)\n", d$snp[which.min(d$p1)], min(d$p1)))
cat(sprintf("OA    min p: %s (%.3g)\n", d$snp[which.min(d$p2)], min(d$p2)))
r <- a$results
cat(sprintf("max SNP.PP.H4: %s (%.3f)\n", r$snp[which.max(r$SNP.PP.H4)],
            max(r$SNP.PP.H4)))
cat(sprintf("SNPs with SNP.PP.H4 > 0.5 : %d\n", sum(r$SNP.PP.H4 > 0.5)))
cat(sprintf("SNPs with SNP.PP.H4 > 0.2 : %d\n", sum(r$SNP.PP.H4 > 0.2)))
