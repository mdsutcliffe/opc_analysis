# Subtype opc90
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(grid)

# source("./functions/build_speciesConversion.R")
source("./functions/collapseIsoforms.R")
source("./functions/normalizeTPM.R")

load("./build/speciesConversion.RData")

f.opc90 <- "./data/rsem_opc90.csv"
f.opc90.info <- "./data/info_opc90.csv"

opc90 <- read.csv(file = f.opc90,stringsAsFactors = F)
opc90.info <- read.csv(file = f.opc90.info,stringsAsFactors = F)

conversion <- convertMouseToHuman(opc90$symbol,human,mouse)

opc90_tpm <- normalizeTPM(rsem = opc90,index_counts = 10:ncol(opc90))

opc90_tpm_human <- opc90_tpm
opc90_tpm_human$symbol <- conversion$HGNC.symbol[match(opc90_tpm_human$symbol,conversion$MGI.symbol)]
opc90_tpm_human <- opc90_tpm_human[complete.cases(opc90_tpm_human),]
opc90_tpm_human <- collapseIsoforms(rsem = opc90_tpm_human,index_counts = 10:ncol(opc90_tpm_human))

opc90_tpm_human_log2 <- cbind(opc90_tpm_human[,1:9],log2(opc90_tpm_human[,10:ncol(opc90_tpm_human)]/100 + 1))

# Load PCs from PMID 31554641
pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv")[,1:10]

opc90_tpm_human_log2_pc <- opc90_tpm_human_log2[match(pc$Gene,opc90_tpm_human_log2$symbol),]

pc1_opc90 <- apply(X = opc90_tpm_human_log2_pc[,10:ncol(opc90_tpm_human_log2_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
pc2_opc90 <- apply(X = opc90_tpm_human_log2_pc[,10:ncol(opc90_tpm_human_log2_pc)],MARGIN = 2,FUN = function(x) {
  pc1 <- pc$PC1 * x
  x <- x - pc1
  pc2 <- pc$PC2 * x
  return(pc2)
})

PC1_opc90 <- colSums(x = pc1_opc90,na.rm = T)
PC2_opc90 <- colSums(x = pc2_opc90,na.rm = T)

row.names(pc1_opc90) <- pc$Gene
row.names(pc2_opc90) <- pc$Gene

pc1_opc90_tencell <- pc1_opc90[complete.cases(pc1_opc90),which(opc90.info$type == "ten-cell")]
pc1_opc90_tencell_annotation <- data.frame(row.names = colnames(pc1_opc90_tencell),
                                           sex = opc90.info$sex[opc90.info$type == "ten-cell"])

pc1_opc90_tencell_top50_inds <- order(rowSums(pc1_opc90_tencell),decreasing = F)[c(1:50,(nrow(pc1_opc90_tencell)-49):nrow(pc1_opc90_tencell))]
pc1_opc90_tencell_top50 <- pc1_opc90_tencell[pc1_opc90_tencell_top50_inds,]

# Plot
{
  # PC 2 vs. PC 1
  pdf(file = "./plots/PC_opc90.pdf",width = 2,height = 2,pointsize = 6)
  par(mar=c(3.5,3.5,1,1))
  plot(x = c(),y = c(),xlim = c(-60,60),ylim = c(0,140),xlab = "PC1",ylab = "PC2",frame = F,axes = F,xaxs = "i",yaxs = "i",mgp = c(2.5,2.5,0))
  axis(side = 1,lwd = 0.5)
  axis(side = 2,lwd = 0.5,las = 1)
  points(x = PC1_opc90[which(opc90.info$type == "ten-cell" & opc90.info$sex == "female")],y = PC2_opc90[which(opc90.info$type == "ten-cell" & opc90.info$sex == "female")],pch = 16,col = "#e41a1c")
  points(x = PC1_opc90[which(opc90.info$type == "ten-cell" & opc90.info$sex == "male")],y = PC2_opc90[which(opc90.info$type == "ten-cell" & opc90.info$sex == "male")],pch = 16,col = "#377eb8")
  legend("topright",legend = c("Female","Male"),pch = c(16,16),col = c("#e41a1c","#377eb8"),box.lwd = 0.5)
  dev.off()
  
  # Histogram of PC 1
  pdf(file = "./plots/PC_histogram_opc90.pdf",width = 2,height = 2,pointsize = 6)
  par(mar=c(3.5,3.5,1,1))
  plot(x = c(),y = c(),xlim = c(-30,30),ylim = c(0,10),xlab = "PC1",ylab = "Number of samples",frame = F,axes = F,xaxs ="i",yaxs = "i",mgp = c(2.5,2.5,0))
  axis(side = 1,at = seq(-30,30,10),lwd = 0.5)
  axis(side = 2,at = seq(0,10,2),lwd = 0.5,las = 1)
  hist(x = PC1_opc90[which(opc90.info$type == "ten-cell")],breaks = seq(-30,30,3),col = "#636363",axes = F,lty = "blank",add = T)
  dev.off()
  
  # Histograms of PC 1 by sex
  pdf(file = "./plots/PC_histogram_opc90_sex.pdf",width = 2,height = 2,pointsize = 6)
  par(mar=c(3.5,3.5,1,1))
  plot(x = c(),y = c(),xlim = c(-30,30),ylim = c(0,10),xlab = "PC1",ylab = "Number of samples",frame = F,axes = F,xaxs ="i",yaxs = "i",mgp = c(2.5,2.5,0))
  axis(side = 1,at = seq(-30,30,10),lwd = 0.5)
  axis(side = 2,at = seq(0,10,2),lwd = 0.5,las = 1)
  replicate(10, { # Hack to show the colors accurately when overlapped
    hist(x = PC1_opc90[which(opc90.info$type == "ten-cell" & opc90.info$sex == "female")],breaks = seq(-30,30,3),col = "#e41a1c40",axes = F,lty = "blank",add = T)
    hist(x = PC1_opc90[which(opc90.info$type == "ten-cell" & opc90.info$sex == "male")],breaks = seq(-30,30,3),col = "#377eb840",axes = F,lty = "blank",add = T)
  })
  legend("topright",legend = c("Female","Male"),pch = c(15,15),col = c("#e41a1c","#377eb8"),box.lwd = 0.5)
  dev.off()
  
  # Heatmap of all PC1 projections
  pdf(file = "./plots/PC_heatmap_opc90.pdf",width = 8,height = 6)
  hm_pc1_opc90 <- pheatmap(mat = t(pc1_opc90_tencell[nrow(pc1_opc90_tencell):1,]),
                           color = rev(brewer.pal(11,"RdBu")),
                           breaks = seq(-3,3,length.out = 12),
                           border_color = NA,
                           cluster_cols = F,
                           clustering_method = "ward.D2",
                           show_colnames = F,
                           show_rownames = F,
                           annotation_row = pc1_opc90_tencell_annotation)
  grid.text(label = "PC1",y = 0.025)
  grid.text(label = "10-cell samples",x = 0.875,rot = 90)
  dev.off()
  
  # Heatmap of top 50 PC1 projections
  pdf(file = "./plots/PC_heatmap_opc90_top50.pdf",width = 8,height = 6)
  hm_pc1_opc90_top50 <- pheatmap(mat = t(pc1_opc90_tencell_top50),
                                 color = rev(brewer.pal(11,"RdBu")),
                                 breaks = seq(-3,3,length.out = 12),
                                 border_color = NA,
                                 cluster_cols = F,
                                 clustering_method = "ward.D2",
                                 show_colnames = F,
                                 show_rownames = F,
                                 annotation_row = pc1_opc90_tencell_annotation,
                                 gaps_col = rep(50,3))
  grid.text(label = "Most influential PC1 genes",y = 0.025)
  grid.text(label = "10-cell samples",x = 0.875,rot = 90)
  dev.off()
  
  # Heatmap of column-scaled top 50 PC1 projections
  pdf(file = "./plots/PC_heatmap_opc90_top50_scale.pdf",width = 8,height = 6)
  hm_pc1_opc90_top50_scale <- pheatmap(mat = t(pc1_opc90_tencell_top50[,hm_pc1_opc90_top50$tree_row$order]),
                                       color = rev(brewer.pal(11,"RdBu")),
                                       breaks = seq(-4,4,length.out = 12),
                                       border_color = NA,
                                       cluster_cols = F,
                                       clustering_method = "ward.D2",
                                       show_colnames = F,
                                       show_rownames = F,
                                       annotation_row = pc1_opc90_tencell_annotation,
                                       gaps_col = rep(50,3),
                                       scale = "column",
                                       clustering_distance_rows = dist(t(pc1_opc90_tencell_top50[,hm_pc1_opc90_top50$tree_row$order])))
  grid.text(label = "Most influential PC1 genes",y = 0.025)
  grid.text(label = "10-cell samples",x = 0.875,rot = 90)
  grid.text(label = "column-scaled",x = 0.925,y = 0.775,rot = 90)
  dev.off()
  
  # Heatmap of absolute values of top 50 PC1 projections
  pdf(file = "./plots/PC_heatmap_opc90_top50_abs.pdf",width = 8,height = 6)
  hm_pc1_opc90_top50 <- pheatmap(mat = abs(t(pc1_opc90_tencell_top50)),
                                 color = brewer.pal(9,"Reds"),
                                 breaks = seq(0,3,length.out = 10),
                                 border_color = NA,
                                 cluster_cols = F,
                                 clustering_method = "ward.D2",
                                 show_colnames = F,
                                 show_rownames = F,
                                 annotation_row = pc1_opc90_tencell_annotation,
                                 gaps_col = rep(50,3))
  grid.text(label = "Most influential PC1 genes",y = 0.025)
  grid.text(label = "10-cell samples",x = 0.875,rot = 90)
  grid.text(label = "absolute value",x = 0.925,y = 0.775,rot = 90)
  dev.off()
  
  # Heatmap of column-scaled absolute values of top 50 PC1 projections
  pdf(file = "./plots/PC_heatmap_opc90_top50_abs_scaled.pdf",width = 8,height = 6)
  hm_pc1_opc90_top50_scale <- pheatmap(mat = abs(t(pc1_opc90_tencell_top50[,hm_pc1_opc90_top50$tree_row$order])),
                                       color = rev(brewer.pal(11,"RdBu")),
                                       breaks = seq(-4,4,length.out = 12),
                                       border_color = NA,
                                       cluster_cols = F,
                                       clustering_method = "ward.D2",
                                       show_colnames = F,
                                       show_rownames = F,
                                       annotation_row = pc1_opc90_tencell_annotation,
                                       gaps_col = rep(50,3),
                                       scale = "column",
                                       clustering_distance_rows = dist(t(pc1_opc90_tencell_top50[,hm_pc1_opc90_top50$tree_row$order])))
  grid.text(label = "Most influential PC1 genes",y = 0.025)
  grid.text(label = "10-cell samples",x = 0.875,rot = 90)
  grid.text(label = "column-scaled\nabsolute value",x = 0.94,y = 0.775,rot = 90)
  dev.off()
}
