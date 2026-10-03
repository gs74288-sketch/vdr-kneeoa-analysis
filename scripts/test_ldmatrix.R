suppressMessages({library(ieugwasr); library(data.table)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
eq <- fread(file.path(base, "data/eqtl_region.tsv"))
eq <- eq[!is.na(rsid) & rsid != "." & !is.na(beta) & !is.na(se) & se>0 & !is.na(p)]
eq <- eq[!duplicated(rsid)]
cand <- eq[p < 5e-5]; cand[, F := (beta/se)^2]; cand <- cand[F > 10]
setorder(cand, p)
top <- head(cand, 200)$rsid
cat("testing ld_matrix on", length(top), "variants\n")
m <- tryCatch(ld_matrix(top, with_alleles=TRUE, pop="EUR"),
              error=function(e){cat("ld_matrix ERR:", substr(conditionMessage(e),1,160), "\n"); NULL})
if (!is.null(m)) {
  cat("ld_matrix OK: dim =", paste(dim(m), collapse=" x "), "\n")
  fwrite(as.data.table(m, keep.rownames="snp"), file.path(base, "outputs/ld_matrix_test.csv"))
  cat("first 3x3:\n"); print(m[1:3,1:3])
}
