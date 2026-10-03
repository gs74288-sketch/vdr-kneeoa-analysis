suppressMessages(library(LDlinkR))
tok <- Sys.getenv("LDLINK_TOKEN")
cat("token present:", nchar(tok) > 0, "| length:", nchar(tok), "\n")
m <- tryCatch(
  LDmatrix(snps = c("rs7975232", "rs7485057", "rs10783222"), pop = "EUR", r2d = "r2", token = tok),
  error = function(e) { cat("LDmatrix ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(m)) { cat("LDmatrix OK\n"); print(m) }
