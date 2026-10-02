suppressMessages({library(ieugwasr); library(data.table)})
cat("1) gwasinfo\n"); print(tryCatch(gwasinfo("eqtl-a-ENSG00000111424"), error=function(e) conditionMessage(e)))
cat("\n2) tophits default (5e-8) on VDR eQTL\n")
t1 <- tryCatch(tophits(id="eqtl-a-ENSG00000111424"), error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(t1)) cat("rows:", nrow(t1), " cols:", paste(names(t1),collapse=","), "\n")
cat("\n3) associations with 11 instruments on VDR eQTL\n")
snps <- c("rs10783222","rs12943633","rs149110519","rs2594836","rs28498283","rs35979828",
          "rs372883","rs4142441","rs58067800","rs7485057","rs7975232")
a <- tryCatch(associations(variants=snps, id="eqtl-a-ENSG00000111424"), error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(a)){ cat("rows:",nrow(a),"\n"); print(head(a,3)) }
cat("\n4) associations same SNPs on knee OA\n")
b <- tryCatch(associations(variants=snps, id="ebi-a-GCST007090"), error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(b)){ cat("rows:",nrow(b),"\n"); print(head(b,3)) }
