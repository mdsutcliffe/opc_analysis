# Subtype tumors
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(pheatmap)
library(RColorBrewer)

# source("./functions/tumorSubtyping.R")
source("./functions/collapseIsoforms.R")
source("./functions/normalizeTPM.R")

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

bulk_collapse$symbol <- conversion$HGNC.symbol[match(bulk_collapse$symbol,conversion$MGI.symbol)]
bulk_collapse <- bulk_collapse[complete.cases(bulk_collapse),]

bulk_collapse <- collapseIsoforms(bulk_collapse,10:ncol(bulk_collapse))

bulk_collapse_tpm <- normalizeTPM(rsem = bulk_collapse,index_counts = 10:ncol(bulk_collapse))

bulk_collapse_tpm_subtype <- bulk_collapse_tpm[,c(3,1,10:ncol(bulk_collapse_tpm))]
names(bulk_collapse_tpm_subtype)[1:2] <- c("NAME","Description")

# write.table(x = "#1.2",file = "./ssgsea.GBM.classification/opcBulk.gct",quote = F,sep = "\t",row.names = F,col.names = F)
# write.table(x = paste(nrow(bulk_collapse_tpm_subtype),ncol(bulk_collapse_tpm_subtype)-2,sep = "\t"),
#             file = "./ssgsea.GBM.classification/opcBulk.gct",quote = F,sep = "\t",append = T,row.names = F,col.names = F)
# write.table(x = bulk_collapse_tpm_subtype,file = "./ssgsea.GBM.classification/opcBulk.gct",quote = F,sep = "\t",row.names = F,append = T)
# 
# 
# source("./ssgsea.GBM.classification/R/msig.library.12.R")
# source("./ssgsea.GBM.classification/R/runSsGSEAwithPermutationR3.R")
# 
# runSsGSEAwithPermutation(profile_data_file = "./ssgsea.GBM.classification/opcBulk.gct",number_perms = 100)
# 
# p_opcBulk <- read.table("./ssgsea.GBM.classification/p_result_opcBulk.gct.txt",header = T, row.names = 1)
# 
# p_opcBulk_tumor <- p_opcBulk[bulk_collapse.info$day == 150,4:6] < 0.05
# 
# phm_annotation <- data.frame(row.names = bulk_collapse.info$name,
#                              genotype = bulk_collapse.info$genotype,
#                              sex = bulk_collapse.info$sex)
# 
# phm_annotation_colors <- list(genotype = c(WT = "#e41a1c",CKO = "#377eb8"),
#                               sex = c(female = "#006d2c",male = "#54278f"))
# 
# pdf(file = "~/subtype_bulk150.pdf",width = 4,height = 6)
# pheatmap(1-p_opcBulk_tumor,color = c("#000000","#FFFFFF"),breaks = c(0,0.05,1),annotation_row = phm_annotation,
#          annotation_colors = phm_annotation_colors,show_rownames = F,main = "150 dpi",cluster_cols = F,cluster_rows = F)
# dev.off()
# 
# p_opcBulk_90 <- p_opcBulk[bulk_collapse.info$day == 90,4:6] < 0.05
# pdf(file = "~/subtype_bulk90.pdf",width = 4,height = 6)
# pheatmap(1-p_opcBulk_90,color = c("#000000","#FFFFFF"),breaks = c(0,0.05,1),annotation_row = phm_annotation,
#          annotation_colors = phm_annotation_colors,show_rownames = F,main = "90 dpi",cluster_cols = F,cluster_rows = F)
# dev.off()
# p_opcBulk_12 <- p_opcBulk[bulk_collapse.info$day == 12,4:6] < 0.05
# 
# pdf(file = "~/subtype_bulk12.pdf",width = 4,height = 6)
# pheatmap(1-p_opcBulk_12,color = c("#000000","#FFFFFF"),breaks = c(0,0.05,1),annotation_row = phm_annotation,
#          annotation_colors = phm_annotation_colors,show_rownames = F,main = "12 dpi",cluster_cols = F,cluster_rows = F)
# dev.off()

conversion_subtype <- convertHumanToMouse(tumorSubtype$GeneSymbol,human,mouse)

bulk_collapse_tpm_subtype_genes <- bulk_collapse_tpm_subtype[match(tumorSubtype$GeneSymbol,bulk_collapse_tpm_subtype$NAME),]
bulk_collapse_tpm_subtype_genes <- bulk_collapse_tpm_subtype_genes[complete.cases(bulk_collapse_tpm_subtype_genes),]
row.names(bulk_collapse_tpm_subtype_genes) <- bulk_collapse_tpm_subtype_genes$NAME
bulk_collapse_tpm_subtype_genes <- bulk_collapse_tpm_subtype_genes[,3:ncol(bulk_collapse_tpm_subtype_genes)]

annotation <- data.frame(row.names = bulk_collapse.info$name,genotype = bulk_collapse.info$genotype)

pdf(file = "plots/subtype_bulk_CKO_scaled.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                       which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                       which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         scale = "row",cluster_rows = F,cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),gaps_col = rep(c(49,93),each=4),gaps_row = c(7,14))
dev.off()

pdf(file = "plots/subtype_bulk_WT_scaled.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"))] + 1)),
         scale = "row",cluster_rows = F,cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),gaps_col = rep(c(49,93),each=4),gaps_row = c(5,9))
dev.off()

pdf(file = "plots/subtype_bulk_CKO.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         cluster_rows = F,cluster_cols = F,color = brewer.pal(9,"Reds"),gaps_col = rep(c(49,93),each=4),gaps_row = c(7,14))
dev.off()

pdf(file = "plots/subtype_bulk_WT.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"))] + 1)),
         cluster_rows = F,cluster_cols = F,color = brewer.pal(9,"Reds"),gaps_col = rep(c(49,93),each=4),gaps_row = c(5,9))
dev.off()

pdf(file = "plots/subtype_bulk.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         cluster_rows = F,cluster_cols = F,color = brewer.pal(9,"Reds"),gaps_col = rep(c(49,93),each=4),gaps_row = c(5,12,16,23,29))
dev.off()

pdf(file = "plots/subtype_bulk_scaled.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         scale = "row",cluster_rows = F,cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),gaps_col = rep(c(49,93),each=4),gaps_row = c(5,12,16,23,29))
dev.off()

pdf(file = "plots/subtype_bulk_tumor.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         cluster_rows = F,cluster_cols = F,show_rownames = F,
         color = brewer.pal(9,"Greys"),gaps_col = rep(c(49,93),each=4),gaps_row = 6,
         annotation_row = annotation,annotation_colors = list(genotype = c(WT = "#31a354",CKO = "#756bb1")))
dev.off()

pdf(file = "plots/subtype_bulk_tumor_scaled.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_genes[,c(which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         scale = "row",cluster_rows = F,show_rownames = F,cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),gaps_col = rep(c(49,93),each=4),gaps_row = 6,
         annotation_row = annotation,annotation_colors = list(genotype = c(WT = "#31a354",CKO = "#756bb1")))
dev.off()

pGSCsig <- c("PDGFRA","OLIG2","SOX10","ASCL1","DLL3","BCAN")
mGSCsig <- c("PTGDS","CD99","FN1","EGFR","VIM","CHI3L1","CD44")
pGSCsig_nuclear <- c("OLIG2","HES6","BCAN","DLL3","OLIG1","SOX10")
mGSCsig_nuclear <- c("TNC","CHI3L1","CD44","CHL1","MEG3","GFAP")

GSC_sigGenes <- c(unique(c(pGSCsig,pGSCsig_nuclear)),unique(c(mGSCsig,mGSCsig_nuclear)))
bulk_collapse_tpm_subtype_GSC <- bulk_collapse_tpm_subtype[match(GSC_sigGenes,bulk_collapse_tpm_subtype$NAME),]
bulk_collapse_tpm_subtype_GSC <- bulk_collapse_tpm_subtype_GSC[complete.cases(bulk_collapse_tpm_subtype_GSC),]
row.names(bulk_collapse_tpm_subtype_GSC) <- bulk_collapse_tpm_subtype_GSC$NAME
bulk_collapse_tpm_subtype_GSC <- bulk_collapse_tpm_subtype_GSC[,3:ncol(bulk_collapse_tpm_subtype_GSC)]

pdf(file = "plots/gscSig_bulk_CKO_scaled.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_GSC[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         scale = "row",cluster_rows = F,cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),gaps_col = rep(length(unique(c(pGSCsig,pGSCsig_nuclear))),each=4),gaps_row = c(7,14))
dev.off()

pdf(file = "plots/gscSig_bulk_CKO.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_GSC[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                       which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                       which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         cluster_rows = F,cluster_cols = F,color = brewer.pal(9,"Greys"),gaps_col = rep(length(unique(c(pGSCsig,pGSCsig_nuclear))),each=4),gaps_row = c(7,14))
dev.off()

pdf(file = "plots/gscSig_bulk.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_GSC[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         cluster_rows = F,cluster_cols = F,color = brewer.pal(9,"Greys"),gaps_col = rep(length(unique(c(pGSCsig,pGSCsig_nuclear))),each=4),gaps_row = c(5,12,16,23,29))
dev.off()

pdf(file = "plots/gscSig_bulk_scaled.pdf",width = 20,height = 12)
pheatmap(mat = t(log2(bulk_collapse_tpm_subtype_GSC[,c(which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 90 & bulk_collapse.info$genotype == "CKO"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"),
                                                         which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "CKO"))] + 1)),
         scale = "row",cluster_rows = F,cluster_cols = F,color = rev(brewer.pal(11,"RdBu")),gaps_col = rep(length(unique(c(pGSCsig,pGSCsig_nuclear))),each=4),gaps_row = c(5,12,16,23,29))
dev.off()