suppressMessages({library(pdftools)})
dir <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14"
pdfs <- list.files(dir, pattern="\\.pdf$", full.names=TRUE)
for (f in pdfs) {
  cat("\n=====", basename(f), "=====\n")
  txt <- tryCatch(pdf_text(f), error=function(e){cat("err\n"); NULL})
  if (is.null(txt)) next
  txt <- paste(txt, collapse=" ")
  gse <- unique(regmatches(txt, gregexpr("GSE[0-9]{4,7}", txt))[[1]])
  if (length(gse)) cat("GSE:", paste(gse, collapse=", "), "\n") else cat("no GSE found\n")
  m <- regmatches(txt, gregexpr("(E-MTAB-[0-9]+|PRJ[A-Z]+[0-9]+|h5ad|GEO|ArrayExpress|Zenodo|figshare)", txt, ignore.case=TRUE))[[1]]
  if (length(m)) cat("other:", paste(unique(m), collapse=", "), "\n")
}
