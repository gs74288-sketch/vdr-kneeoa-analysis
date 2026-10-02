suppressMessages({library(Seurat); library(data.table); library(Matrix)})
out <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/outputs"
merged <- readRDS(file.path(out, "GSE152805_seurat.rds"))
merged <- JoinLayers(merged)
cat("cells:", ncol(merged), " cluster levels:", length(levels(merged$seurat_clusters)), "\n")

mk <- list(HomC=c("ACAN","COL2A1","SOX9","COMP","PRG4","TIMP3"),
           HTC=c("COL10A1","IBSP","MMP13","RUNX2","IHH","ALPL"),
           FC=c("COL1A1","COL3A1","DCN","FN1"),
           ProC=c("MKI67","TOP2A","UBE2C"),
           EC=c("VEGFA","HIF1A","ANGPTL4"))
data <- LayerData(merged, assay="RNA", layer="data")
genes <- unique(unlist(mk)); genes <- genes[genes %in% rownames(data)]

cl <- as.character(merged$seurat_clusters)
avg <- sapply(sort(unique(cl), decreasing=FALSE), function(c) rowMeans(data[genes, cl==c, drop=FALSE]))
avg <- as.data.table(as.data.frame(avg), keep.rownames="gene")
fwrite(avg, file.path(out, "GSE152805_cluster_markers.csv"))

clusters <- colnames(avg)[-1]
M <- as.matrix(avg[, -1]); rownames(M) <- avg$gene
# z-score each gene ACROSS clusters, then average per marker set
Z <- t(scale(t(M)))
Z[is.na(Z)] <- 0
csum <- data.table(cluster=clusters)
for (k in names(mk)) {
  g <- intersect(mk[[k]], rownames(Z))
  csum[[k]] <- if (length(g)) colMeans(Z[g, , drop=FALSE]) else 0
}
csum[, state := names(mk)[max.col(as.matrix(.SD), ties.method="first")], .SDcols=names(mk)]
print(csum)
fwrite(csum, file.path(out, "GSE152805_cluster_state_assignment.csv"))

stmap <- setNames(csum$state, csum$cluster)
md <- data.table(cell=colnames(merged), cluster=cl, donor=merged$donor, site=merged$site,
                 nCount=merged$nCount_RNA)
md[, state := stmap[cluster]]
md[, vdr := data["VDR", ]]
md[, detected := as.integer(vdr > 0)]
md[, total := merged$nCount_RNA]

st <- md[, .(n=.N, detected=sum(detected), frac=mean(detected)), by=state][order(-frac)]
print(st); fwrite(st, file.path(out, "GSE152805_VDR_by_state.csv"))

reg <- c("HomC")
tab <- table(factor(md$state %in% reg, c(TRUE,FALSE)), md$detected)
ft <- fisher.test(tab)
cat("\nHomC(-like) vs other states: Fisher OR =", round(ft$estimate,3),
    " P =", format.pval(ft$p.value), "\n")
pd <- md[, .(n=.N, detected=sum(detected), frac=mean(detected)),
         by=.(donor, state)][order(donor, -frac)]
fwrite(pd, file.path(out, "GSE152805_VDR_per_donor_state.csv")); print(dcast(pd, donor ~ state, value.var="frac"))

md[, is_reg := as.integer(state %in% reg)]
m <- glm(detected ~ is_reg + log1p(total) + donor, family=quasibinomial(), data=md)
cat("\n--- quasibinomial: detected ~ state + log(depth) + donor ---\n")
print(summary(m)$coefficients)
# OR for is_reg
b <- summary(m)$coefficients["is_reg", "Estimate"]; s <- summary(m)$coefficients["is_reg","Std. Error"]
cat("OR for HomC(-like):", round(exp(b),3), "95% CI", round(exp(b-1.96*s),3), "-", round(exp(b+1.96*s),3),
    " p =", format.pval(summary(m)$coefficients["is_reg","Pr(>|t|)"]), "\n")
cat("\nDONE_GSE152805_ANALYSIS\n")
