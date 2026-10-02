suppressMessages({library(ieugwasr); library(data.table)})
jwt <- Sys.getenv("OPENGWAS_JWT")
cat("env OPENGWAS_JWT present:", nchar(jwt) > 0, " len:", nchar(jwt), "\n")
r <- tryCatch(gwasinfo("ebi-a-GCST007090"), error=function(e){cat("gwasinfo ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(r)) { cat("gwasinfo OK ->", nrow(r), "row(s)\n"); print(as.data.frame(r)[, intersect(c("id","trait","sample_size","ncase","ncontrol","population"), names(r))]) }
