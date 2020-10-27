#Figure 3B

library(dplyr)
library('edgeR')
library(pheatmap)
library(RColorBrewer)
library(rstudioapi)

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

#Start with the raw read counts of ALL samples aligned using Hisat
rawdata<-read.table("gene_read_counts_table_all_final.tsv",header=TRUE, stringsAsFactors=FALSE, row.names=1)
mapping<-read.table("ENSG_ID2Name.txt", header=FALSE, stringsAsFactors=FALSE, row.names=1)

#Check dimensions
dim(rawdata)

#Append genenames to rawdata and reorder to make 1st column
rawdata$symbol<-mapping$V2
rawdata<- rawdata[c(25,1:24)]

#There are duplicate gene names in all these outputs, this function sums them to prevent issues downstream:
sum_duplicates <- function(tpmdata)
{
  index <- (2:ncol(tpmdata))
  duplicated_genes <- tpmdata[duplicated(tpmdata$symbol),] #all the repeats
  unique_genes <- unique(duplicated_genes$symbol) #just the name of the repeated genes ONCE
  # Sum duplicate gene names
  for (iGene in unique_genes) {
    toSum <- tpmdata[which(tpmdata$symbol == iGene),]
    indexOfMaxCount <- which.max(rowSums(toSum[,index]))
    summed <- colSums(toSum[,index])
    tpmdata[names(indexOfMaxCount),index] <- summed
    toDelete <- toSum[-indexOfMaxCount,]
    
    tpmdata <- tpmdata[!rownames(tpmdata) %in% rownames(toDelete),]
  }
  return(tpmdata)
}
rawdata_nodups <- sum_duplicates(rawdata)

#Get genes that are expressed at greater than 5 TPM in all samples 
#(analyzed in another script but saved gene names for future use)

filtered_genes<-read.table("Genes_expressed>5TPMallsamples.txt",header=TRUE)

filt_counts<-rawdata_nodups[rawdata_nodups$symbol %in% filtered_genes$x,]

#EdgeR analysis
gene_names<-filt_counts$symbol

#Perform THREE differential gene expression analyses because edgeR can only do differential
#gene expression b/w 2 classes 
filt_5E=filt_counts[,c(3,5,7,9,10:13)] #Get shGFP 2019 samples and shNRF2v1
#filt_5E=filt_counts[,c(3,5,7,9,14:17)] #Get shGFP 2019 samples and shNRF2v4
#filt_DCIS=filt_counts[,c(18:25)] #Get DCIS shCon and shNRF2v1 samples
class<-factor(c(rep("shCon",4),rep("shNRF2",4)))
genes=rownames(filt_DCIS)

y <- DGEList(counts=filt_DCIS, genes=genes, group=class)
nrow(y)

# TMM Normalization
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# Differential expression test
et <- exactTest(y)

# Print top genes
topTags(et)

# Print number of up/down significant genes at FDR = 0.05  significance level
summary(de <- decideTestsDGE(et, p=.05))
detags <- rownames(y)[as.logical(de)]

# Output DE genes
# Matrix of significantly DE genes
mat <- cbind(
  genes,gene_names,
  sprintf('%0.3f',log10(et$table$PValue)),
  sprintf('%0.3f',et$table$logFC)
)[as.logical(de),]
colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")

# Order by log fold change
o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
mat <- mat[o,]

#Save table
#write.table(mat, file="DE_genes_5EshGFPvsshNRF2v1-9-25-19.txt", quote=FALSE, row.names=FALSE, sep="\t")
#write.table(mat, file="DE_genes_5EshGFPvsshNRF2v4-9-25-19.txt", quote=FALSE, row.names=FALSE, sep="\t")
#write.table(mat, file="DE_genes_DCISshLacZv3vsshNRF2v1-9-25-19.txt", quote=FALSE, row.names=FALSE, sep="\t")

#Analyze DE genes ------------
tab_5Ev1=read.table(file="DE_genes_5EshGFPvsshNRF2v1-9-25-19.txt",header=TRUE, row.names=NULL, sep="\t")
tab_5Ev4=read.table(file="DE_genes_5EshGFPvsshNRF2v4-9-25-19.txt",header=TRUE, row.names=NULL, sep="\t")
tab_DCIS=read.table(file="DE_genes_DCISshLacZv3vsshNRF2v1-9-25-19.txt",header=TRUE, row.names=NULL, sep="\t")

#Find DE genes that change in the same direction between 5E shGFP vs shNRF2v1 and 5E shGFP vs shNRF2v4-
sigup_5Ev1<-subset(tab_5Ev1,Log_fold_change>0)
sigdown_5Ev1<-subset(tab_5Ev1,Log_fold_change<0)

sigup_5Ev4<-subset(tab_5Ev4,Log_fold_change>0)
sigdown_5Ev4<-subset(tab_5Ev4,Log_fold_change<0)

sigup_5E<-intersect(sigup_5Ev1$Gene_Name,sigup_5Ev4$Gene_Name)
sigdown_5E<-intersect(sigdown_5Ev1$Gene_Name,sigdown_5Ev4$Gene_Name)

sigup_DCIS<-subset(tab_DCIS,Log_fold_change>0)
sigdown_DCIS<-subset(tab_DCIS,Log_fold_change<0)

#Find genes in common b/w 5E and DCIS
common_up<-data.frame(intersect(sigup_5E,sigup_DCIS$Gene_Name))
common_down<-data.frame(intersect(sigdown_5E,sigdown_DCIS$Gene_Name))

#Get gene names so I can do MSigDb analysis with them
#write.table(common_up, file="Upreg_common_genes_5EshNRF2v1v4andDCISshNRF2v1.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, sep="\t")
#write.table(common_down, file="Common_downreg_genes_5EshNRF2v1v4andDCISshNRF2v1.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, sep="\t")

#MSigDb analysis- getting DE genenames that are in the BRCA1,CHEK2,ATM networks 
#for annotation in figure 3B
library(cmapR)

BRCA1<-parse.grp('BRCA1_PCC.grp')
BRCA1overlap<-data.frame(intersect(BRCA1,common_up$intersect.sigup_5E..sigup_DCIS.Gene_Name.))
colnames(BRCA1overlap)<-"OverlapBRCA1_commonup"
#write.table(BRCA1overlap, file="Upreg_common_genes_inBRCA1_PCCnetwork.txt",quote=FALSE, col.names=TRUE,row.names=FALSE, sep="\t")

CHEK2<-parse.grp('CHEK2_PCC.grp')
CHEK2overlap<-data.frame(intersect(CHEK2,common_up$intersect.sigup_5E..sigup_DCIS.Gene_Name.))
colnames(CHEK2overlap)<-"OverlapCHEK2_commonup"
#write.table(CHEK2overlap, file="Upreg_common_genes_inCHEK2_PCCnetwork.txt",quote=FALSE, col.names=TRUE,row.names=FALSE, sep="\t")

ATM<-parse.grp('ATM_PCC.grp')
ATMoverlap<-data.frame(intersect(ATM,common_up$intersect.sigup_5E..sigup_DCIS.Gene_Name.))
colnames(ATMoverlap)<-"OverlapATM_commonup"
#write.table(ATMoverlap, file="Upreg_common_genes_inATM_PCCnetwork.txt",quote=FALSE, col.names=TRUE,row.names=FALSE, sep="\t")

#Find unique genes
BRCA1unique<-setdiff(BRCA1overlap$OverlapBRCA1_commonup,CHEK2overlap$OverlapCHEK2_commonup)
BRCA1unique<-setdiff(BRCA1unique,ATMoverlap$OverlapATM_commonup)

CHEK2unique<-setdiff(CHEK2overlap$OverlapCHEK2_commonup,BRCA1overlap$OverlapBRCA1_commonup)



#Pull out DE genes (commonup+commondown) from TPM dataset to make heatmap
colnames(common_down)<-"gene"
colnames(common_up)<-"gene"
DEgenes<-rbind(common_down,common_up)

#write.table(DEgenes, file="Common_genes_5EshNRF2v1v4andDCISshNRF2v1.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, sep="\t")

#Use TPM data for generating heatmap
liz_data <- read.csv('5EshGFPshNRF2v1v4DCIS.comshLacZv3shNRF2v1_TPM_geneid.csv',stringsAsFactors = F)

#There are duplicate gene names in all these outputs, this function sums them to prevent issues downstream:
sum_duplicates <- function(tpmdata)
{
  index <- (3:ncol(tpmdata))
  duplicated_genes <- tpmdata[duplicated(tpmdata$symbol),] #all the repeats
  unique_genes <- unique(duplicated_genes$symbol) #just the name of the repeated genes ONCE
  # Sum duplicate gene names
  for (iGene in unique_genes) {
    toSum <- tpmdata[which(tpmdata$symbol == iGene),]
    indexOfMaxCount <- which.max(rowSums(toSum[,index]))
    summed <- colSums(toSum[,index])
    tpmdata[names(indexOfMaxCount),index] <- summed
    toDelete <- toSum[-indexOfMaxCount,]
    
    tpmdata <- tpmdata[!rownames(tpmdata) %in% rownames(toDelete),]
  }
  return(tpmdata)
}
liz_data <- sum_duplicates(liz_data)

#For ease of downstream steps, at this point i would set the gene names as the row names and maybe only get data columns

rownames(liz_data) <- liz_data$symbol
liz_data_mat <- liz_data[,3:22]

#Try >5TPM in all samples
threshold <- 5
proportion<-1
liz_data_filt_allsamples<-liz_data_mat[rowSums(liz_data_mat[,]>threshold) == (proportion*ncol(liz_data_mat)),]

liz_data_filt<-liz_data_filt_allsamples

#Log2 transform
liz_data_filt  <- log2(liz_data_filt +1)

#Pull out DEgenes
DEgenesmat<-liz_data_filt[match(DEgenes$gene,row.names(liz_data_filt)),]
#write.csv(DEgenesmat,'Matrix_commonupreganddownreg_5EDCISshNRF2-9-25-19.csv',row.names = T)

mybreaks <- seq(-4, 4, length.out=8)

p<-pheatmap(DEgenesmat , scale='row', clustering_distance_rows = "euclidean",cluster_cols =F,
            breaks=mybreaks,clustering_method = "ward.D2", show_rownames = F, labels_row = F,
            col=rev(brewer.pal(7 ,"RdBu")))