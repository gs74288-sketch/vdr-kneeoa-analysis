suppressMessages({library(Seurat); library(data.table)})
out <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/outputs"
merged <- readRDS(file.path(out, "GSE152805_seurat.rds"))
merged <- JoinLayers(merged)
data <- LayerData(merged, assay="RNA", layer="data")
cl <- as.character(merged$seurat_clusters)
clusters <- sort(unique(cl))
cat("cells:", ncol(merged), "clusters:", length(clusters), "\n")

panel <- list(
  HomC  = c("ACAN","COL2A1","SOX9","COMP","PRG4","TIMP3","COL11A1"),
  HTC   = c("COL10A1","IBSP","MMP13","RUNX2","IHH","ALPL","PHOSPHO1"),
  preHTC= c("COL10A1","COL9A1","SPP1"),
  FC    = c("COL1A1","COL3A1","DCN","FN1","LUM"),
  ProC  = c("MKI67","TOP2A","UBE2C","CENPF","CCNB1"),
  EC    = c("VEGFA","HIF1A","ANGPTL4","SLC2A1"),
  RegC  = c("CCL3","CCL4","CCL4L2","CXCL8","CXCL1","CXCL2","CXCL3","CXCL6","IL6","NFKBIA",
            "TNFAIP3","SOCS3","CEBPB","FOS","FOSB","JUN","JUNB","ATF3","NFKB1","RELA","ICAM1","CXCL5"))
genes <- unique(unlist(panel)); genes <- genes[genes %in% rownames(data)]
cat("panel genes found:", length(genes), "/", length(unique(unlist(panel))), "\n")
cat("VDR present:", "VDR" %in% rownames(data), "\n")

mean_by_cl <- sapply(clusters, function(c) rowMeans(data[genes, cl==c, drop=FALSE]))
colnames(mean_by_cl) <- clusters
Z <- t(scale(t(mean_by_cl))); Z[is.na(Z)] <- 0

sc <- data.table(cluster=clusters, n=as.integer(table(cl)[clusters]))
for (k in names(panel)) {
  g <- intersect(panel[[k]], rownames(Z))
  sc[[k]] <- if (length(g)) colMeans(Z[g,,drop=FALSE]) else 0
}
# VDR detection per cluster
det <- sapply(clusters, function(c) mean(data["VDR", cl==c] > 0))
sc[, vdr_frac := det[cluster]]
sc[, vdr_n := as.integer(sapply(clusters, function(c) sum(data["VDR", cl==c] > 0)))]
fwrite(sc, file.path(out, "GSE152805_cluster_scores_regC.csv"))
cat("\nPer-cluster scores & VDR detection (sorted by RegC score):\n")
print(sc[order(-RegC), .(cluster, n, HomC=round(HomC,2), HTC=round(HTC,2), FC=round(FC,2),
                         ProC=round(ProC,2), EC=round(EC,2), RegC=round(RegC,2),
                         vdr_n, vdr_frac=round(vdr_frac,4))], nrows=25)

# define RegC-like = clusters in top 25% of RegC score and with RegC score > 0.5
thr <- quantile(sc$RegC, 0.75)
regc_like <- sc[RegC >= max(0.5, thr)]$cluster
cat("\nRegC-like clusters:", paste(regc_like, collapse=","), "\n")

md <- data.table(cell=colnames(merged), cluster=cl, d=as.integer(data["VDR",] > 0),
                 total=merged$nCount_RNA, donor=merged$donor)
md[, regc := cluster %in% regc_like]
cat("\nVDR detection in RegC-like vs other:\n")
print(md[, .(n=.N, detected=sum(d), frac=round(mean(d),4)), by=regc])
ft <- fisher.test(table(md$regc, md$d))
cat("Fisher OR =", round(ft$estimate,3), " P =", format.pval(ft$p.value), "\n")
m <- glm(d ~ regc + log1p(total) + donor, family=quasibinomial(), data=md)
cf <- summary(m)$coefficients["regcTRUE",]
cat("Adjusted OR =", round(exp(cf[1]),3), "95% CI", round(exp(cf[1]-1.96*cf[2]),3), "-",
    round(exp(cf[1]+1.96*cf[2]),3), " p =", format.pval(cf[4]), "\n")
# correlation of cluster-level RegC score vs VDR detection
cat("\nSpearman(cluster RegC score, cluster VDR frac) =",
    round(cor(sc$RegC, sc$vdr_frac, method="spearman"),3), "\n")
cat("Spearman(cluster HomC score, cluster VDR frac) =",
    round(cor(sc$HomC, sc$vdr_frac, method="spearman"),3), "\n")
cat("\nDONE_REGC\n")
