# Whole transcriptome clustering to identify outliers

library(pheatmap)
library(RColorBrewer)
library(ggbiplot)

source("./functions/normalizeTPM.R")

# File paths
bulk.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/rsem_opcBulk.csv"
bulk.info.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/info_opcBulk.csv"

# Import
bulk <- read.csv(bulk.path,stringsAsFactors = F)
bulk.info <- read.csv(bulk.info.path,stringsAsFactors = F)

# Assign row names
row.names(bulk.info) <- bulk.info$name

# Keep only tdT+ samples
bulk <- bulk[,c(1:9,9+which(bulk.info$celltype == "tdTpositive"))]
bulk.info <- bulk.info[bulk.info$celltype == "tdTpositive",]

# TPM normalize and Log2 transform data
bulk_tpm <- normalizeTPM(bulk,10:ncol(bulk))
row.names(bulk_tpm) <- bulk_tpm$symbol
bulk_tpm <- bulk_tpm[,10:ncol(bulk_tpm)]
bulk_tpm_log <- log2(bulk_tpm + 1)

# Separate by dpi
bulk_tpm_log_12 <- bulk_tpm_log[,bulk.info$day == 12]
bulk_tpm_log_90 <- bulk_tpm_log[,bulk.info$day == 90]
bulk_tpm_log_150 <- bulk_tpm_log[,bulk.info$day == 150]

# Pull out genes where half the samples > 0
bulk_tpm_log_12_filter <- bulk_tpm_log_12[apply(bulk_tpm_log_12,1,function(x) sum(x > 0) > floor(length(x)/2)),]
bulk_tpm_log_90_filter <- bulk_tpm_log_90[apply(bulk_tpm_log_90,1,function(x) sum(x > 0) > floor(length(x)/2)),]
bulk_tpm_log_150_filter <- bulk_tpm_log_150[apply(bulk_tpm_log_150,1,function(x) sum(x > 0) > floor(length(x)/2)),]


# Clustering + heatmaps
plotHeatmap <- function(rsem_tpm_log, info) {
  annotation <- info[,c("genotype","sex")]
  annotation_colors <- list(genotype = c(WT = "#BDBDBD",CKO = "#636363"),sex = c(female = "#41AB5D",male = "#807DBA"))
  phm <- pheatmap(mat = rsem_tpm_log,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  breaks = seq(from = -4,to = 4,length.out = 12),
                  scale = "row",
                  show_rownames = F,
                  color = rev(brewer.pal(11,"RdBu")),
                  annotation_col = annotation,
                  annotation_colors = annotation_colors)
  return(phm)
}

phm_12 <- plotHeatmap(rsem_tpm_log = bulk_tpm_log_12_filter,
                      info = bulk.info[bulk.info$day == 12,])

phm_90 <- plotHeatmap(rsem_tpm_log = bulk_tpm_log_90_filter,
                      info = bulk.info[bulk.info$day == 90,])
outliers <- c("opcBulk_90dpi_CKO_20647_tdTpositive_GAGCT","opcBulk_90dpi_WT_21167_tdTpositive")
phm_90_noOutliers <- plotHeatmap(rsem_tpm_log = bulk_tpm_log_90_filter[,!(names(bulk_tpm_log_90_filter) %in% outliers)],
                                 info = bulk.info[bulk.info$day == 90,])

phm_150 <- plotHeatmap(rsem_tpm_log = bulk_tpm_log_150_filter,
                       info = bulk.info[bulk.info$day == 150,])

# PCA
pca_12 <- prcomp(x = t(bulk_tpm_log_12_filter),center = T,scale. = T)
pca_plot_12 <- ggbiplot(pcobj = pca_12,
                      labels = names(bulk_tpm_log_12_filter),
                      var.axes = F,
                      groups = bulk.info$sex[bulk.info$day == 12]) +
  scale_color_manual(values=c("#de2d26","#3182bd"))
pca_plot_12

pca_90 <- prcomp(x = t(bulk_tpm_log_90_filter),center = T,scale. = T)
pca_plot_90 <- ggbiplot(pcobj = pca_90,
                        labels = names(bulk_tpm_log_90_filter),
                        var.axes = F,
                        groups = bulk.info$sex[bulk.info$day == 90]) +
  scale_color_manual(values=c("#de2d26","#3182bd"))
pca_plot_90

pca_90_noOutliers <- prcomp(x = t(bulk_tpm_log_90_filter[,!(names(bulk_tpm_log_90_filter) %in% outliers)]),center = T,scale. = T)
pca_plot_90_noOutliers <- ggbiplot(pcobj = pca_90_noOutliers,
                        labels = names(bulk_tpm_log_90_filter[,!(names(bulk_tpm_log_90_filter) %in% outliers)]),
                        var.axes = F,
                        groups = bulk.info$sex[bulk.info$day == 90 & !(bulk.info$name %in% outliers)]) +
  scale_color_manual(values=c("#de2d26","#3182bd"))
pca_plot_90_noOutliers

pca_150 <- prcomp(x = t(bulk_tpm_log_150_filter),center = T,scale. = T)
pca_plot_150 <- ggbiplot(pcobj = pca_150,
                        labels = names(bulk_tpm_log_150_filter),
                        var.axes = F,
                        groups = bulk.info$sex[bulk.info$day == 150]) +
  scale_color_manual(values=c("#de2d26","#3182bd"))
pca_plot_150
