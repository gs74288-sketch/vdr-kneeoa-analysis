suppressMessages({library(ieugwasr); library(data.table)})
jwt <- Sys.getenv("OPENGWAS_JWT")
cat("JWT present:", nchar(jwt) > 0, " length:", nchar(jwt), "\n")
r <- tryCatch(gwasinfo("ebi-a-GCST007090"), error=function(e){cat("gwasinfo ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(r)) { cat("gwasinfo OK\n"); print(r) }
eq <- tryCatch(tophits(id="eqtl-a-ENSG00000111424", pval=1), error=function(e){cat("tophits ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(eq)) { cat("VDR cis eQTL SNPs:", nrow(eq), "\n"); cat("cols:", paste(names(eq), collapse=","), "\n"); print(head(eq,3)) }
