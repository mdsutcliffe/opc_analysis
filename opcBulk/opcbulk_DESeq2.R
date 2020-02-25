# DESeq2 on bulk RNA-sequencing data

library(pheatmap)
library(RColorBrewer)
library(ggbiplot)
library(DESeq2)
library(hypeR)

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

# Convert all columns to factors
bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],factor)

# Re-level to WT
bulk_collapse.info$genotype <- relevel(x = bulk_collapse.info$genotype,ref = "WT")

# TPM normalize and Log2 transform data for later analysis
bulk_collapse_tpm <- normalizeTPM(rsem = bulk_collapse,index_counts = 10:ncol(bulk_collapse))
row.names(bulk_collapse_tpm) <- bulk_collapse_tpm$symbol
bulk_collapse_tpm <- bulk_collapse_tpm[,10:ncol(bulk_collapse_tpm)]
bulk_collapse_tpm_log <- log2(bulk_collapse_tpm + 1)

# Get counts matrix and remove gene annotation
row.names(bulk_collapse) <- bulk_collapse$symbol
bulk_collapse <- bulk_collapse[,10:ncol(bulk_collapse)]

# Convert all columns to factors
bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],factor)

# For DESeq2 compatability
row.names(bulk_collapse.info) <- bulk_collapse.info$name

# DESeq2 requires integer counts
bulk_collapse <- round(bulk_collapse)

# DESeq2
runDESeq2 <- function(rsem,info,day,design) {
  rsem <- rsem[,info$day == day]
  info <- info[info$day == day,]
  
  rsem <- as.matrix(rsem)
  mode(rsem) <- "integer"
  
  modelMatrix <- model.matrix(object = design,data = info)
  
  rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
  
  dds <- DESeqDataSetFromMatrix(countData = rsem_filter,colData = info,design = modelMatrix)
  dds <- DESeq(object = dds,quiet = T,test = "LRT",reduced = model.matrix(object = formula(~ sex),data = info))
  
  res <- results(object = dds,contrast = list("genotypeCKO"))
  resOrdered <- res[order(res$padj),]
  
  genesDE <- row.names(resOrdered)[!is.na(resOrdered$padj) & resOrdered$padj < 0.05]
  
  message(paste("design =", format(design),"\n  ",sprintf(fmt = "%6s",day),"dpi :",sprintf(fmt = "%4s",length(genesDE)),"genes"))
  
  return(list("design" = design,
              "modelMatrix" = modelMatrix,
              "genesDE" = genesDE,
              "results" = resOrdered,
              "dds" = dds))
}

design <- formula(~ genotype * sex)
res_12 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 12,design = design)
res_90 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 90,design = design)
res_150 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 150,design = design)

runDESeq2 <- function(rsem,info,day,design) {
  rsem <- rsem[,info$day == day]
  info <- info[info$day == day,]
  
  rsem <- as.matrix(rsem)
  mode(rsem) <- "integer"
  
  modelMatrix <- model.matrix(object = design,data = info)
  
  rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
  
  dds <- DESeqDataSetFromMatrix(countData = rsem_filter,colData = info,design = modelMatrix)
  dds <- DESeq(object = dds,quiet = T,test = "LRT",reduced = model.matrix(object = formula(~ sex),data = info))
  
  res <- results(object = dds,contrast = list("genotypeCKO"))
  resOrdered <- res[order(res$padj),]
  
  genesDE <- row.names(resOrdered)[!is.na(resOrdered$padj) & resOrdered$padj < 0.05]
  
  message(paste("design =", format(design),"\n  ",sprintf(fmt = "%6s",day),"dpi :",sprintf(fmt = "%4s",length(genesDE)),"genes"))
  
  return(list("design" = design,
              "modelMatrix" = modelMatrix,
              "genesDE" = genesDE,
              "results" = resOrdered,
              "dds" = dds))
}

design <- formula(~ genotype * sex)
res_12 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 12,design = design)
res_90 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 90,design = design)
res_150 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 150,design = design)

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
                  show_colnames = F,
                  color = rev(brewer.pal(11,"RdBu")),
                  annotation_col = annotation,
                  annotation_colors = annotation_colors,
                  main = paste(nrow(rsem_tpm_log),"DE genes"))
  return(phm)
}


# plotHeatmap(bulk_collapse_tpm_log[res_12$genesDE,bulk_collapse.info$day == 12],bulk_collapse.info)
# plotHeatmap(bulk_collapse_tpm_log[res_90$genesDE,bulk_collapse.info$day == 90],bulk_collapse.info)
# plotHeatmap(bulk_collapse_tpm_log[res_150$genesDE,bulk_collapse.info$day == 150],bulk_collapse.info)


design <- formula(~ genotype + sex)
res_plus_12 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 12,design = design)
res_plus_90 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 90,design = design)
res_plus_150 <- runDESeq2(rsem = bulk_collapse,info = bulk_collapse.info,day = 150,design = design)

# plotHeatmap(bulk_collapse_tpm_log[res_plus_12$genesDE,bulk_collapse.info$day == 12],bulk_collapse.info)
# plotHeatmap(bulk_collapse_tpm_log[res_plus_90$genesDE,bulk_collapse.info$day == 90],bulk_collapse.info)
# plotHeatmap(bulk_collapse_tpm_log[res_plus_150$genesDE,bulk_collapse.info$day == 150],bulk_collapse.info)


# Enrichment

msigdb_path <- msigdb_download_all(species = "Mus musculus",output_dir = "./external")
hallmark <- msigdb_fetch(msigdb_path = msigdb_path,symbol = "H")

bulk_150_up <- hypeR(signature = row.names(res_plus_150$results[res_plus_150$results$log2FoldChange > 0 & !is.na(res_plus_150$results$padj) & res_plus_150$results$padj < 0.05,]),gsets = hallmark,fdr_cutoff = 0.01,do_plots = T)
pdf("./plots/enrichment_bulk_150_up.pdf",width = 8,height = 5)
hyp_dots(bulk_150_up,top = 100)
dev.off()

bulk_150_down <- hypeR(signature = row.names(res_plus_150$results[res_plus_150$results$log2FoldChange < 0 & !is.na(res_plus_150$results$padj) & res_plus_150$results$padj < 0.05,]),gsets = hallmark,fdr_cutoff = 0.05,do_plots = T)
pdf("./plots/enrichment_bulk_150_down.pdf",width = 8,height = 5)
hyp_dots(bulk_150_down,top = 100)
dev.off()


geneTable150 <- res_plus_150$results[,c("log2FoldChange","pvalue","padj")]
write.csv(x = geneTable150,file = "./temp/bulk150_genes.csv",quote = F)

geneTable12 <- res_plus_12$results[,c("log2FoldChange","pvalue","padj")]
write.csv(x = geneTable12,file = "./temp/bulk12_genes.csv",quote = F)

geneTable90 <- res_plus_90$results[,c("log2FoldChange","pvalue","padj")]
write.csv(x = geneTable90,file = "./temp/bulk90_genes.csv",quote = F)


# Volcano plot
fc <- res_plus_150$results$log2FoldChange[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj >= 0.05]
p <- res_plus_150$results$pvalue[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj >= 0.05]

fc_up <- res_plus_150$results$log2FoldChange[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange > 0]
p_up <- res_plus_150$results$pvalue[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange > 0]

fc_down <- res_plus_150$results$log2FoldChange[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange < 0]
p_down <- res_plus_150$results$pvalue[!is.na(res_plus_150$results$pvalue) & res_plus_150$results$padj < 0.05 & res_plus_150$results$log2FoldChange < 0]

png("./plots/volcano_150.png",width = 2000,height = 2000,res = 300)
par(mar=c(4.1,4.1,0.5,0.5))
plot(fc,-log(p,base = 10),
     xlim = c(-10,10),
     ylim = c(0,100),
     xlab = expression("Log"[2]*" fold change"),
     ylab = expression("-Log"[10]*"("*italic("p")*"-value)"),
     las = 1,frame = F)
points(fc_up,-log(p_up,base = 10),col = "#de2d26")
points(fc_down,-log(p_down,base = 10),col = "#3182bd")
legend(x = "topright",legend = c("Increased","Decreased"),col = c("#de2d26","#3182bd"),pch=c(1,1))
dev.off()

plotCounts(dds = res_plus_150$dds,gene = "Pbx3",intgroup = "genotype")
plotCounts(dds = res_plus_12$dds,gene = "Rpl30",intgroup = "genotype")

wt <- bulk_collapse_tpm_log["Rpl30",which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT")]
cko <- bulk_collapse_tpm_log["Rpl30",which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO")]
png(filename = "./plots/plotCounts_bulk12_Rpl30.png",width = 800,height = 1200,res = 250)
par(mar=c(2.1,4.1,2.5,0.5))
plot(c(rep(1,length(wt)),rep(2,length(cko))),c(wt,cko),
     xlim = c(0.8,2.2),ylim = c(0,12),las = 1,xlab = "",axes = F,ylab = expression("Log"[2]*"(TPM + 1)"),main = "12 dpi bulk\nRpl30")
axis(1,at = c(1,2),labels = c("WT","CKO"))
axis(2,las = 1)
dev.off()

wt <- bulk_collapse_tpm_log["Rpl34",which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT")]
cko <- bulk_collapse_tpm_log["Rpl34",which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO")]
png(filename = "./plots/plotCounts_bulk12_Rpl34.png",width = 800,height = 1200,res = 250)
par(mar=c(2.1,4.1,2.5,0.5))
plot(c(rep(1,length(wt)),rep(2,length(cko))),c(wt,cko),
     xlim = c(0.8,2.2),ylim = c(0,12),las = 1,xlab = "",axes = F,ylab = expression("Log"[2]*"(TPM + 1)"),main = "12 dpi bulk\nRpl34")
axis(1,at = c(1,2),labels = c("WT","CKO"))
axis(2,las = 1)
dev.off()

wt <- bulk_collapse_tpm_log["Rplp1",which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT")]
cko <- bulk_collapse_tpm_log["Rplp1",which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO")]
png(filename = "./plots/plotCounts_bulk12_Rplp1.png",width = 800,height = 1200,res = 250)
par(mar=c(2.1,4.1,2.5,0.5))
plot(c(rep(1,length(wt)),rep(2,length(cko))),c(wt,cko),
     xlim = c(0.8,2.2),ylim = c(0,12),las = 1,xlab = "",axes = F,ylab = expression("Log"[2]*"(TPM + 1)"),main = "12 dpi bulk\nRplp1")
axis(1,at = c(1,2),labels = c("WT","CKO"))
axis(2,las = 1)
dev.off()

# See if pGSC/mGSC genes exist within DE genes

load("./temp/speciesConversion.RData")

# File paths
bulk.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/rsem_opcBulk.csv"
bulk.info.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/info_opcBulk.csv"

# Import
bulk <- read.csv(bulk.path,stringsAsFactors = F)
bulk.info <- read.csv(bulk.info.path,stringsAsFactors = F)

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

# Convert all columns to factors
bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],factor)

conversion <- convertMouseToHuman(bulk_collapse$symbol,human,mouse)


# bulk_collapse$symbol <- conversion$HGNC.symbol[match(bulk_collapse$symbol,conversion$MGI.symbol)]
# bulk_collapse <- bulk_collapse[complete.cases(bulk_collapse),]

bulk_collapse <- collapseIsoforms(bulk_collapse,10:ncol(bulk_collapse))

bulk_collapse_tpm <- normalizeTPM(rsem = bulk_collapse,index_counts = 10:ncol(bulk_collapse))

pc <- read.csv("./temp/219026_2_supp_5782669_py1bdv.csv")[,1:3]
pc_convert <- convertHumanToMouse(pc$Gene,human,mouse)
pc_mouse <- unique(pc_convert$MGI.symbol)

intersect(res_plus_150$genesDE,pc_mouse)
df <- data.frame(x = c(rep("DE",length(res_plus_150$genesDE)),rep("PC",length(pc_mouse))),
                 y = c(res_plus_150$genesDE,pc_mouse))
overlap.test(df,10175,1000)

#2C

bulk_2c <- bulk_collapse_tpm_log[,which(bulk_collapse.info$genotype == "WT" & bulk_collapse.info$day == 150)]
opc_markers <- c("Pdgfra","Lnx1","Dcn","Mmp15","Cdo1")
pericyte_markers <- c("Fmod","Rps2","Igf2","Gpc3","Ogn")
microglia_markers <- c("Slfn2","Gpr84","Ccr7","Bcl2a1d","Tnf")
neuron_markers <- c("Reln","Nhlh2","Slc17a6","Trp73","Lhx5")
astrocyte_markers <- c("Hgf","Aqp4","Itih3","Bmpr1b","Itga7")
endothelial_markers <- c("Cldn5","Ttr","Ly6a","Madcam1","Akr1c14")
nfo_markers <- c("Gp1bb","Tmem108","Fyn","Ust","Mical3")
mo_markers <- c("Gjb1","Ndrg1","Ppp1r14a","Adssl1","Aspa")

bulk_2c[opc_markers,]
bulk_2c[pericyte_markers,]
bulk_2c[microglia_markers,]
bulk_2c[neuron_markers,]
bulk_2c[astrocyte_markers,]

bulk_2c <- bulk_2c[c(opc_markers,pericyte_markers,nfo_markers,mo_markers,microglia_markers,neuron_markers,astrocyte_markers,endothelial_markers),]
pheatmap(t(bulk_2c),cluster_cols = F)

png(filename = "./plots/figure_2c.png",width = 1000,height = 500,res = 100)
boxplot(t(bulk_2c),at = (1:55)[(1:55) %% 7 != 0 & (1:55) %% 7 != 6],ylim = c(0,10),xlim = c(0,55),
        main = "Bulk 150 dpi wild-type",ylab = expression("Log"[2]*"(TPM + 1)"),axes=F,frame = F)
axis(2,las = 1)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5,55.5),labels = F)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5)+3.25,lwd=0,labels = c("OPC","NFO","MO","pericyte","microglia","neuron","astrocyte","endothelial"))
dev.off()


bulk_2c <- bulk_collapse_tpm_log[,which(bulk_collapse.info$genotype == "WT" & bulk_collapse.info$day == 12)]
bulk_2c <- bulk_2c[c(opc_markers,pericyte_markers,nfo_markers,mo_markers,microglia_markers,neuron_markers,astrocyte_markers,endothelial_markers),]
png(filename = "./plots/figure_2c_12dpiWT.png",width = 1000,height = 500,res = 100)
boxplot(t(bulk_2c),at = (1:55)[(1:55) %% 7 != 0 & (1:55) %% 7 != 6],ylim = c(0,10),xlim = c(0,55),
        main = "Bulk 12 dpi wild-type",ylab = expression("Log"[2]*"(TPM + 1)"),axes=F,frame = F)
axis(2,las = 1)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5,55.5),labels = F)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5)+3.25,lwd=0,labels = c("OPC","NFO","MO","pericyte","microglia","neuron","astrocyte","endothelial"))
dev.off()

geneOrder <- unlist(lapply(list(opc_markers,pericyte_markers,nfo_markers,mo_markers,microglia_markers,neuron_markers,astrocyte_markers,endothelial_markers),function(x) {
  sort_inds <- order(apply(bulk_2c[x,],1,median),decreasing = T)
  x[sort_inds]
}))

bulk_2c <- bulk_collapse_tpm_log[,which(bulk_collapse.info$genotype == "WT" & bulk_collapse.info$day == 150)]

geneOrder <- unlist(lapply(list(opc_markers,pericyte_markers,nfo_markers,mo_markers,microglia_markers,neuron_markers,astrocyte_markers,endothelial_markers),function(x) {
  sort_inds <- order(apply(bulk_2c[x,],1,median),decreasing = T)
  x[sort_inds]
}))

bulk_2c <- bulk_2c[geneOrder,]
png(filename = "./plots/figure_2c_ordered.png",width = 1000,height = 500,res = 100)
boxplot(t(bulk_2c),at = (1:55)[(1:55) %% 7 != 0 & (1:55) %% 7 != 6],ylim = c(0,10),xlim = c(0,55),
        main = "Bulk 150 dpi wild-type",ylab = expression("Log"[2]*"(TPM + 1)"),axes=F,frame = F,col = rep(c("#aaaaaaaa",rep("#FFFFFFFF",3),"#aaaaaaaa"),8),border = rep(c("#aaaaaa",rep("#000000",3),"#AAAAAA"),8))
axis(2,las = 1)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5,55.5),labels = F)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5)+3.25,lwd=0,labels = c("OPC","NFO","MO","pericyte","microglia","neuron","astrocyte","endothelial"))
dev.off()

bulk_2c <- bulk_collapse_tpm_log[,which(bulk_collapse.info$genotype == "WT" & bulk_collapse.info$day == 12)]
geneOrder <- unlist(lapply(list(opc_markers,pericyte_markers,nfo_markers,mo_markers,microglia_markers,neuron_markers,astrocyte_markers,endothelial_markers),function(x) {
  sort_inds <- order(apply(bulk_2c[x,],1,median),decreasing = T)
  x[sort_inds]
}))
bulk_2c <- bulk_2c[geneOrder,]
png(filename = "./plots/figure_2c_ordered_12dpiWT.png",width = 1000,height = 500,res = 100)
boxplot(t(bulk_2c),at = (1:55)[(1:55) %% 7 != 0 & (1:55) %% 7 != 6],ylim = c(0,10),xlim = c(0,55),
        main = "Bulk 12 dpi wild-type",ylab = expression("Log"[2]*"(TPM + 1)"),axes=F,frame = F,col = rep(c("#aaaaaaaa",rep("#FFFFFFFF",3),"#aaaaaaaa"),8),border = rep(c("#aaaaaa",rep("#000000",3),"#AAAAAA"),8))
axis(2,las = 1)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5,55.5),labels = F)
axis(1,at = c(0,6.5,13.5,20.5,27.5,34.5,41.5,48.5)+3.25,lwd=0,labels = c("OPC","NFO","MO","pericyte","microglia","neuron","astrocyte","endothelial"))
dev.off()













