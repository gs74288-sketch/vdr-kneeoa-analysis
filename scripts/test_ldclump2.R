suppressMessages({library(ieugwasr); library(data.table)})
base <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002"
eq <- fread(file.path(base, "data/eqtl_region.tsv"))
eq <- eq[!is.na(rsid) & rsid != "." & !is.na(beta) & !is.na(se) & se>0 & !is.na(p)]
eq <- eq[!duplicated(rsid)]
cand <- eq[p < 5e-5]; cand[, F := (beta/se)^2]; cand <- cand[F > 10]
cat("candidates:", nrow(cand), "\n")
for (i in 1:3) {
  cl <- tryCatch(ld_clump(data.frame(rsid=cand$rsid, pval=cand$p, id="eqtl-a-ENSG00000111424"),
                          clump_kb=10000, clump_r2=0.001, pop="EUR"),
                 error=function(e){cat("attempt", i, "ERR:", substr(conditionMessage(e),1,120), "\n"); NULL})
  if (!is.null(cl)) {
    cat("SUCCESS: LD-clumped instruments =", nrow(cl), "\n")
    print(cl)
    fwrite(as.data.table(cl), file.path(base, "outputs/cis_ldclumped_instruments.csv"))
    break
  }
  Sys.sleep(5)
}
