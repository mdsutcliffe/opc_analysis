# Whole transcriptome clustering - outliers removed and replicates collapsed

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

# Remove outliers
outliers <- c("opcBulk_90dpi_CKO_20647_tdTpositive_GAGCT","opcBulk_90dpi_WT_21167_tdTpositive")
bulk <- bulk[,!(names(bulk) %in% outliers)]
bulk.info <- bulk.info[!(bulk.info$name %in% outliers),]

# Rename samples if they're from the same mouse as an outlier
for (iOutlier in outliers) {
  if (length(strsplit(x = iOutlier,split = "_")[[1]]) == 6) {
    iRootName <- paste(strsplit(x = iOutlier,split = "_")[[1]][-6],collapse = "_")
    names(bulk)[which(grepl(iRootName,names(bulk)))] <- iRootName
    bulk.info$name[which(grepl(iRootName,bulk.info$name))] <- iRootName
  }
}

bulk_collapse <- bulk
bulk_collapse.info <- bulk.info

# Convert all columns to characters
bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],as.character)

# Collapse replicates from same mouse
bulk_collapse_unique_replicateIndex.info <- unique(bulk_collapse.info[,c("day","genotype","sex","brainID")])
indsToRemove <- c()
for (i in 1:nrow(bulk_collapse_unique_replicateIndex.info)) {
  inds <- which(bulk_collapse.info$day == bulk_collapse_unique_replicateIndex.info$day[i] &
                  bulk_collapse.info$genotype == bulk_collapse_unique_replicateIndex.info$genotype[i] &
                  bulk_collapse.info$sex == bulk_collapse_unique_replicateIndex.info$sex[i] &
                  bulk_collapse.info$brainID == bulk_collapse_unique_replicateIndex.info$brainID[i])
  if (length(inds) > 1) {
    bulk_collapse[,9+inds[1]] <- rowSums(bulk_collapse[,9+inds])
    bulk_collapse.info[inds[1],"nCells"] <- paste(bulk_collapse.info[inds,"nCells"],collapse="|")
    bulk_collapse.info[inds[1],"runID"] <- paste(bulk_collapse.info[inds,"runID"],collapse="|")
    bulk_collapse.info[inds[1],"barcode"] <- paste(bulk_collapse.info[inds,"barcode"],collapse="|")
    bulk_collapse.info[inds[1],"replicate"] <- paste(bulk_collapse.info[inds,"replicate"],collapse="|")
    indsToRemove <- c(indsToRemove,inds[2:length(inds)])
    
    iNewName <- paste(strsplit(x = names(bulk_collapse)[9+inds[1]],split = "_")[[1]][-6],collapse = "_")
    names(bulk_collapse)[9+inds[1]] <- iNewName
    bulk_collapse.info$name[inds[1]] <- iNewName
  }
}
bulk_collapse <- bulk_collapse[,-(9+indsToRemove)]
bulk_collapse.info <- bulk_collapse.info[-indsToRemove,]

# TPM normalize and Log2 transform data
bulk_collapse_tpm <- normalizeTPM(bulk_collapse,10:ncol(bulk_collapse))
row.names(bulk_collapse_tpm) <- bulk_collapse_tpm$symbol
bulk_collapse_tpm <- bulk_collapse_tpm[,10:ncol(bulk_collapse_tpm)]
bulk_collapse_tpm_log <- log2(bulk_collapse_tpm + 1)

# Separate by dpi
bulk_collapse_tpm_log_12 <- bulk_collapse_tpm_log[,bulk_collapse.info$day == 12]
bulk_collapse_tpm_log_90 <- bulk_collapse_tpm_log[,bulk_collapse.info$day == 90]
bulk_collapse_tpm_log_150 <- bulk_collapse_tpm_log[,bulk_collapse.info$day == 150]

# Pull out genes where half the samples > 0
bulk_collapse_tpm_log_12_filter <- bulk_collapse_tpm_log_12[apply(bulk_collapse_tpm_log_12,1,function(x) sum(x > 0) > floor(length(x)/2)),]
bulk_collapse_tpm_log_90_filter <- bulk_collapse_tpm_log_90[apply(bulk_collapse_tpm_log_90,1,function(x) sum(x > 0) > floor(length(x)/2)),]
bulk_collapse_tpm_log_150_filter <- bulk_collapse_tpm_log_150[apply(bulk_collapse_tpm_log_150,1,function(x) sum(x > 0) > floor(length(x)/2)),]

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

phm_12 <- plotHeatmap(rsem_tpm_log = bulk_collapse_tpm_log_12_filter,
                      info = bulk_collapse.info[bulk_collapse.info$day == 12,])

phm_90 <- plotHeatmap(rsem_tpm_log = bulk_collapse_tpm_log_90_filter,
                      info = bulk_collapse.info[bulk_collapse.info$day == 90,])

phm_150 <- plotHeatmap(rsem_tpm_log = bulk_collapse_tpm_log_150_filter,
                       info = bulk_collapse.info[bulk_collapse.info$day == 150,])

# PCA
pca_12 <- prcomp(x = t(bulk_collapse_tpm_log_12_filter),center = T,scale. = T)
pca_plot_12 <- ggbiplot(pcobj = pca_12,
                        labels = names(bulk_collapse_tpm_log_12_filter),
                        var.axes = F,
                        groups = bulk_collapse.info$sex[bulk_collapse.info$day == 12]) +
  scale_color_manual(values=c("#de2d26","#3182bd"))
pca_plot_12

pca_90 <- prcomp(x = t(bulk_collapse_tpm_log_90_filter),center = T,scale. = T)
pca_plot_90 <- ggbiplot(pcobj = pca_90,
                        labels = names(bulk_collapse_tpm_log_90_filter),
                        var.axes = F,
                        groups = bulk_collapse.info$sex[bulk_collapse.info$day == 90]) +
  scale_color_manual(values=c("#de2d26","#3182bd"))
pca_plot_90

pca_150 <- prcomp(x = t(bulk_collapse_tpm_log_150_filter),center = T,scale. = T)
pca_plot_150 <- ggbiplot(pcobj = pca_150,
                         labels = names(bulk_collapse_tpm_log_150_filter),
                         var.axes = F,
                         groups = bulk_collapse.info$sex[bulk_collapse.info$day == 150]) +
  scale_color_manual(values=c("#de2d26","#3182bd"))
pca_plot_150