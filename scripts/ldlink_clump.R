suppressMessages({library(data.table); library(LDlinkR)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
out  <- file.path(base, "outputs")
tok  <- Sys.getenv("LDLINK_TOKEN")
stopifnot(nchar(tok) > 0)

eq <- fread(file.path(base, "data/eqtl_region.tsv"))
eq <- eq[!is.na(rsid) & rsid != "." & !is.na(beta) & !is.na(se) & se > 0 & !is.na(p) & !is.na(position)]
eq <- eq[!duplicated(rsid)]
cand <- eq[p < 5e-5]
cand[, F := (beta/se)^2]
cand <- cand[F > 10]
setorder(cand, position)
cat("candidates (P<5e-5, F>10):", nrow(cand), "\n")

ld_mat <- function(snps) {
  for (k in 1:4) {
    m <- tryCatch(LDmatrix(snps = snps, pop = "EUR", r2d = "r2", token = tok),
                  error = function(e) { cat("  LDmatrix err:", substr(conditionMessage(e),1,90), "\n"); NULL })
    if (!is.null(m)) return(m)
    Sys.sleep(8)
  }
  NULL
}

CH <- 250
survivors <- cand[0]
chunks <- split(seq_len(nrow(cand)), ceiling(seq_len(nrow(cand)) / CH))
for (ci in seq_along(chunks)) {
  idx <- chunks[[ci]]
  sub <- cand[idx]
  m <- ld_mat(sub$rsid)
  if (is.null(m)) { cat("chunk", ci, "failed -> keeping all\n"); survivors <- rbind(survivors, sub); next }
  rn <- m[[1]]; M <- as.matrix(m[, -1]); rownames(M) <- rn; colnames(M) <- rn
  keep <- character(0)
  for (s in sub[order(p)]$rsid) {
    if (!(s %in% rownames(M))) next
    if (length(keep) == 0 || all(M[s, keep] < 0.001)) keep <- c(keep, s)
  }
  survivors <- rbind(survivors, sub[rsid %in% keep])
  cat("chunk", ci, "/", length(chunks), ": ", length(keep), "/", nrow(sub), " kept | total", nrow(survivors), "\n")
  Sys.sleep(2)
}

cat("\nsurvivors after chunked clumping:", nrow(survivors), "\n")
if (nrow(survivors) > 1 && nrow(survivors) <= 1000) {
  m <- ld_mat(survivors$rsid)
  if (!is.null(m)) {
    rn <- m[[1]]; M <- as.matrix(m[, -1]); rownames(M) <- rn; colnames(M) <- rn
    keep <- character(0)
    for (s in survivors[order(p)]$rsid) {
      if (!(s %in% rownames(M))) next
      if (length(keep) == 0 || all(M[s, keep] < 0.001)) keep <- c(keep, s)
    }
    survivors <- survivors[rsid %in% keep]
    cat("global pass ->", nrow(survivors), "instruments\n")
  }
}
fwrite(survivors[, .(SNP = rsid, chr, position, ea, nea, eaf, beta, se, p, F)],
       file.path(out, "cis_ldclumped_instruments.csv"))
print(survivors[, .(rsid, position, p, F)])
cat("\nDONE_LDLINK_CLUMP\n")
