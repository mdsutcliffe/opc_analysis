setwd("~/Github/opc_analysis")

library(infercnv)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./new/inferCNV.tpm.matrix",
                                    annotations_file="./new/inferCNV.annotation.txt",
                                    delim="\t",
                                    gene_order_file="./new/inferCNV.gene_order.txt",
                                    ref_group_names=c("bulk12WT","bulk90WT","bulk150WT"))

infercnv_obj = infercnv::run(infercnv_obj,out_dir = "~/Downloads",
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,debug = F)

save.image(file = "~/Github/opc_analysis/new/inferCNV_bulkreference.RData")

{
  sc <- read.csv(file = "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/RNA-seq/Singh et al 2018/janes_gse75330_rsem_counts.csv",stringsAsFactors = F)
  
  source("./functions/normalizeTPM.R")
  
  sc <- normalizeTPM(rsem = sc,index_counts = 10:ncol(sc))
  
  source("./import/import_opc12.R")
  source("./import/import_opc90.R")
  source("./import/import_opcBulk.R")
  
  opc12$tpm <- opc12$tpm[nchar(opc12$tpm$chr) <= 2,]
  opc90$tpm <- opc90$tpm[nchar(opc90$tpm$chr) <= 2,]
  sc <- sc[nchar(sc$chr) <= 2,]
  
  sc <- sc[order(match(sc$chr,c(as.character(1:19),"X","Y","MT")),sc$start),]
  opc12$tpm <- opc12$tpm[order(match(opc12$tpm$chr,c(as.character(1:19),"X","Y","MT")),opc12$tpm$start),]
  opc90$tpm <- opc90$tpm[order(match(opc90$tpm$chr,c(as.character(1:19),"X","Y","MT")),opc90$tpm$start),]
  bulk$tpm <- bulk$tpm[order(match(bulk$tpm$chr,c(as.character(1:19),"X","Y","MT")),bulk$tpm$start),]
  
  commonGenes <- intersect(sc$ensgene,opc12$tpm$ensgene)
  
  tpm <- cbind(data.frame(row.names = opc12$tpm$symbol[match(commonGenes,opc12$tpm$ensgene)]),
               opc12$tpm[match(commonGenes,opc12$tpm$ensgene),9 + which(opc12$info$type == "ten-cell")],
               opc90$tpm[match(commonGenes,opc90$tpm$ensgene),9 + which(opc90$info$type == "ten-cell")],
               bulk$tpm[match(commonGenes,bulk$tpm$ensgene),10:ncol(bulk$tpm)],
               sc[match(commonGenes,sc$ensgene),10:ncol(sc)])
  
  annotation <- cbind(names(tpm),
                      c(rep(x = "opc12",times = 56),
                        rep(x = "opc90",times = 56),
                        paste0("bulk",bulk$info$day,bulk$info$genotype),
                        rep(x = "SC",times = 112)))
  
  gene_order <- sc[match(commonGenes,sc$ensgene),c("symbol","chr","start","end")]
}


library(pheatmap)
library(RColorBrewer)
unique(annotation[,2])
nrow(annotation)
phm_ann <- data.frame(row.names = annotation[,1],group = annotation[,2])
head(phm_ann)
range(as.numeric(infercnv_obj@expr.data))
which(diff(as.numeric(phm_ann$group)) != 0)
png(filename = "./new/cnvplot_bulkreference.png",width = 2000,height = 1500,res = 200)
pheatmap(mat = t(infercnv_obj@expr.data[,!(annotation[,2] %in% c("opc12","opc90","SC"))]),cluster_cols = F,cluster_rows = F,show_colnames = F,show_rownames = F,color = brewer.pal(n = 3,name = "RdBu"),breaks = c(0,0.667,1.5,10),
         annotation_row = phm_ann[annotation[,2] != "SC",1,drop = F],annotation_col = infercnv_obj@gene_order[,1,drop = F],
         gaps_row = rep(x = which(diff(as.numeric(phm_ann$group[!(annotation[,2] %in% c("opc12","opc90","SC"))])) != 0),each = 2))
dev.off()