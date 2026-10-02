suppressMessages({library(Seurat); library(data.table); library(Matrix); library(ggplot2)})
dir <- "C:/Users/liuha/OpenClawVDR/gse152805"
out <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/outputs"

samples <- data.frame(
  gsm = c("GSM4626766","GSM4626767","GSM4626768","GSM4626769","GSM4626770","GSM4626771"),
  tag = c("oLT_113","oLT_116","oLT_118","MT_113","MT_116","MT_118"),
  donor = c("113","116","118","113","116","118"),
  site = c("oLT","oLT","oLT","MT","MT","MT"), stringsAsFactors = FALSE)

objs <- list()
for (i in seq_len(nrow(samples))) {
  s <- samples[i, ]
  mat <- ReadMtx(mtx=file.path(dir, paste0(s$gsm, "_OA_", s$tag, ".matrix.mtx.gz")),
                 cells=file.path(dir, paste0(s$gsm, "_OA_", s$tag, ".barcodes.tsv.gz")),
                 features=file.path(dir, paste0(s$gsm, "_OA_", s$tag, ".genes.tsv.gz")),
                 feature.column=2)
  o <- CreateSeuratObject(mat, project=s$tag, min.cells=3, min.features=200)
  o$sample <- s$tag; o$donor <- s$donor; o$site <- s$site
  objs[[s$tag]] <- o
  cat("loaded", s$tag, "cells:", ncol(o), "\n"); flush.console()
}
merged <- merge(objs[[1]], y=objs[-1], add.cell.ids=samples$tag)
cat("total cells:", ncol(merged), "\n")

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern="^MT-")
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
cat("cells after QC:", ncol(merged), "\n"); flush.console()

merged <- NormalizeData(merged, verbose=FALSE)
merged <- FindVariableFeatures(merged, nfeatures=2000, verbose=FALSE)
merged <- ScaleData(merged, verbose=FALSE)
merged <- RunPCA(merged, npcs=30, verbose=FALSE)
merged <- FindNeighbors(merged, dims=1:20, verbose=FALSE)
merged <- FindClusters(merged, resolution=0.6, verbose=FALSE)
merged <- RunUMAP(merged, dims=1:20, verbose=FALSE)
saveRDS(merged, file.path(out, "GSE152805_seurat.rds"))

# marker-based annotation
mk <- list(HomC=c("ACAN","COL2A1","SOX9","COMP","PRG4","TIMP3"),
           HTC=c("COL10A1","IBSP","MMP13","RUNX2","IHH","ALPL"),
           FC=c("COL1A1","COL3A1","DCN","FN1"),
           ProC=c("MKI67","TOP2A","UBE2C"),
           EC=c("VEGFA","HIF1A","ANGPTL4"))
avg <- AverageExpression(merged, features=unique(unlist(mk)), group.by="seurat_clusters", assays="RNA")$RNA
avg <- as.data.table(as.data.frame(avg), keep.rownames="gene")
fwrite(avg, file.path(out, "GSE152805_cluster_markers.csv"))

score <- function(g) { g <- g[g %in% rownames(merged)]; if(!length(g)) return(rep(0, ncol(merged))); colMeans(GetAssayData(merged, slot="data")[g, , drop=FALSE]) }
for (k in names(mk)) merged[[paste0("score_", k)]] <- score(mk[[k]])
csum <- data.table(cluster=levels(merged$seurat_clusters))
for (k in names(mk)) csum[[k]] <- sapply(csum$cluster, function(cl) mean(merged[[paste0("score_",k)]][merged$seurat_clusters==cl, ]))
print(csum)

# assign each cluster to nearest state by max score
csum[, state := names(mk)[max.col(.SD, ties.method="first")], .SDcols=names(mk)]
ann <- setNames(csum$state, csum$cluster)
merged$state <- ann[as.character(merged$seurat_clusters)]
cat("cluster->state:\n"); print(csum[, .(cluster, state)]); flush.console()
fwrite(csum, file.path(out, "GSE152805_cluster_state_assignment.csv"))

# VDR detection by state
md <- as.data.table(merged@meta.data, keep.rownames="cell")
md[, vdr := GetAssayData(merged, slot="data")["VDR", ]]
md[, detected := as.integer(vdr > 0)]
md[, total := nCount_RNA]
st <- md[, .(n=.N, detected=sum(detected), frac=mean(detected)), by=state][order(-frac)]
print(st); fwrite(st, file.path(out, "GSE152805_VDR_by_state.csv"))

regstates <- c("HomC")
tab <- table(factor(md$state %in% regstates, c(TRUE,FALSE)), md$detected)
ft <- fisher.test(tab)
cat("HomC vs others Fisher OR =", round(ft$estimate,3), " P =", format.pval(ft$p.value), "\n")

per_donor <- md[, .(n=.N, frac=mean(detected)), by=.(donor)][order(donor)]
print(per_donor); fwrite(per_donor, file.path(out, "GSE152805_VDR_per_donor.csv"))

# depth-adjusted quasibinomial with donor covariate
md[, is_reg := as.integer(state %in% regstates)]
m <- glm(detected ~ is_reg + log1p(total) + donor, family=quasibinomial(), data=md)
cat("\n--- depth/donor-adjusted model (VDR ~ state + depth + donor) ---\n"); print(summary(m)$coefficients)

cat("\nDONE_GSE152805\n")
