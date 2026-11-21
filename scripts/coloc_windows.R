## Windowed colocalisation sensitivity: does the answer depend on the interval?
suppressMessages(library(coloc))

IN  <- "C:/Users/liuha/OpenClawVDR/repo/tables/coloc_harmonised_input_GCST90566800.csv"
OUT <- "C:/Users/liuha/OpenClawVDR/repo/tables"
d   <- read.csv(IN, stringsAsFactors = FALSE)

N_EQTL <- 31684
N_CASE <- median(d$ncase, na.rm = TRUE)
N_CTRL <- median(d$ncontrol, na.rm = TRUE)
N_OA   <- N_CASE + N_CTRL
S_OA   <- N_CASE / N_OA

EQTL_PEAK <- 47845054   # rs7975232
OA_PEAK   <- 48009982   # rs11168351

fit <- function(dd, tag) {
  if (nrow(dd) < 50) { cat(sprintf("%-34s too few SNPs (%d)\n", tag, nrow(dd))); return(NULL) }
  d1 <- list(snp = dd$snp, beta = dd$b1, varbeta = dd$se1^2, N = N_EQTL,
             type = "quant", sdY = 1, MAF = dd$maf)
  d2 <- list(snp = dd$snp, beta = dd$b2h, varbeta = dd$se2^2, N = N_OA,
             type = "cc", s = S_OA, MAF = dd$maf)
  r <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5))$summary
  cat(sprintf("%-34s n=%5d  H3=%.4f  H4=%.3e  nGW=%d\n",
              tag, r[["nsnps"]], r[["PP.H3.abf"]], r[["PP.H4.abf"]],
              sum(dd$p2 < 5e-8)))
  data.frame(window = tag, nsnps = r[["nsnps"]], PP.H3 = r[["PP.H3.abf"]],
             PP.H4 = r[["PP.H4.abf"]])
}

rows <- list()
rows[[length(rows) + 1]] <- fit(d, "full region 46.8-48.9 Mb")

for (w in c(10000, 25000, 50000, 100000)) {
  rows[[length(rows) + 1]] <- fit(d[abs(d$pos - EQTL_PEAK) <= w, ],
                                  sprintf("eQTL peak +/- %d kb", w / 1000))
}
for (w in c(25000, 50000, 100000)) {
  rows[[length(rows) + 1]] <- fit(d[abs(d$pos - OA_PEAK) <= w, ],
                                  sprintf("OA peak +/- %d kb", w / 1000))
}
# the two peaks together, narrow
rows[[length(rows) + 1]] <- fit(d[d$pos >= EQTL_PEAK - 25000 & d$pos <= OA_PEAK + 25000, ],
                                "both peaks, 47.82-48.04 Mb")

res <- do.call(rbind, Filter(Negate(is.null), rows))
write.csv(res, file.path(OUT, "coloc_GCST90566800_windows.csv"), row.names = FALSE)
cat("\nwritten coloc_GCST90566800_windows.csv\n")
