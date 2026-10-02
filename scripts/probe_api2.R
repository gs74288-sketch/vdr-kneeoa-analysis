suppressMessages({library(ieugwasr); library(data.table)})
snps <- c("rs10783222","rs12943633","rs149110519","rs2594836","rs28498283","rs35979828",
          "rs372883","rs4142441","rs58067800","rs7485057","rs7975232")
cat("A) associations 11 SNPs on VDR eQTL\n")
a <- tryCatch(associations(variants=snps, id="eqtl-a-ENSG00000111424"), error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(a)){ cat("rows:",nrow(a)," cols:",paste(names(a),collapse=","),"\n"); print(head(a,3)) }
cat("\nB) associations 11 SNPs on knee OA\n")
b <- tryCatch(associations(variants=snps, id="ebi-a-GCST007090"), error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
if(!is.null(b)){ cat("rows:",nrow(b),"\n"); print(head(b,3)) }
