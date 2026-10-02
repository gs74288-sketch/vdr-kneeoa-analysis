suppressMessages({library(ieugwasr); library(data.table)})
ids <- c("finn-b-M13_KNEEARTHROSE","finn-b-M13_KNEEARTHROSE_EXMORE","finn-b-M13_KNEEARTHROSIS",
         "finn-b-M13_ARTHROSIS","finn-b-M13_ARTHROSIS_EXMORE","finn-b-M13_GONARTHROSE",
         "finn-b-M13_PATELLAR_ARTHROSIS","ieu-b-5090")
for (id in ids) {
  r <- tryCatch(gwasinfo(id), error=function(e) NULL)
  if (!is.null(r) && nrow(r) > 0) {
    cat(sprintf("%-38s OK | trait=%s | n=%s | nsnp=%s | cases=%s\n",
                id, r$trait[1], r$sample_size[1], r$nsnp[1],
                ifelse("ncase" %in% names(r), r$ncase[1], "-")))
  } else cat(sprintf("%-38s not found\n", id))
}
