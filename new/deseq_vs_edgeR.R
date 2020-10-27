setwd(dir = "~/Github/opc_analysis/")

library(edgeR)

source("./import/import_opcBulk.R")

# Get counts matrix and remove gene annotation
row.names(bulk$rsem) <- bulk$rsem$symbol
bulk$rsem <- bulk$rsem[,10:ncol(bulk$rsem)]

# Convert all columns to factors
bulk$info[,names(bulk$info)] <- lapply(bulk$info[,names(bulk$info)],factor)

bulk$info$genotype <- relevel(x = bulk$info$genotype,ref = "WT")

# For DESeq2 compatability
row.names(bulk$info) <- bulk$info$name

# DESeq2 requires integer counts
bulk$rsem <- round(bulk$rsem)

# 150 dpi
rsem <- bulk$rsem[,bulk$info$day == 150]
info <- bulk$info[bulk$info$day == 150,c("genotype","sex")]
rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
y <- edgeR::DGEList(counts = rsem_filter,samples = info,group = info$genotype)
design <- model.matrix(object = formula(~ sex + genotype),data = info)
y <- edgeR::estimateCommonDisp(y = y,design = design)
fit <- edgeR::glmFit(y = y,design = design)
lrt <- edgeR::glmLRT(glmfit = fit,coef = 3)
de150 <- edgeR::decideTestsDGE(object = lrt)

# 12 dpi
rsem <- bulk$rsem[,bulk$info$day == 12]
info <- bulk$info[bulk$info$day == 12,c("genotype","sex")]
rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
y <- edgeR::DGEList(counts = rsem_filter,samples = info,group = info$genotype)
design <- model.matrix(object = formula(~ sex + genotype),data = info)
y <- edgeR::estimateCommonDisp(y = y,design = design)
fit <- edgeR::glmFit(y = y,design = design)
lrt <- edgeR::glmLRT(glmfit = fit,coef = 3)
de12 <- edgeR::decideTestsDGE(object = lrt)

# 90 dpi
rsem <- bulk$rsem[,bulk$info$day == 90]
info <- bulk$info[bulk$info$day == 90,c("genotype","sex")]
rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
y <- edgeR::DGEList(counts = rsem_filter,samples = info,group = info$genotype)
design <- model.matrix(object = formula(~ sex + genotype),data = info)
y <- edgeR::estimateCommonDisp(y = y,design = design)
fit <- edgeR::glmFit(y = y,design = design)
lrt <- edgeR::glmLRT(glmfit = fit,coef = 3)
de90 <- edgeR::decideTestsDGE(object = lrt)

# Compare to DESeq results

load("./build/opcBulk_DESeq2_results.RData")

summary(de12)
row.names(de12)[de12 != 0]
length(intersect(bulk$deseq2$de12$genesDE,row.names(de12)[de12 != 0]))
length(setdiff(bulk$deseq2$de12$genesDE,row.names(de12)[de12 != 0]))
length(setdiff(row.names(de12)[de12 != 0],bulk$deseq2$de12$genesDE))

summary(de90)
row.names(de90)[de90 != 0]
length(intersect(bulk$deseq2$de90$genesDE,row.names(de90)[de90 != 0]))
length(setdiff(bulk$deseq2$de90$genesDE,row.names(de90)[de90 != 0]))
length(setdiff(row.names(de90)[de90 != 0],bulk$deseq2$de90$genesDE))

summary(de150)
row.names(de150)[de150 != 0]
length(intersect(bulk$deseq2$de150$genesDE,row.names(de150)[de150 != 0]))
length(setdiff(bulk$deseq2$de150$genesDE,row.names(de150)[de150 != 0]))
length(setdiff(row.names(de150)[de150 != 0],bulk$deseq2$de150$genesDE))

# Plot a heatmap to see if I'm crazy
phm <- bulk$log2[bulk$log2$symbol %in% row.names(de90)[de90 != 0],9+which(bulk$info$day == 90)]
library(pheatmap)
library(RColorBrewer)
head(phm)
pheatmap(mat = phm,scale = "row",color = rev(brewer.pal(n = 11,name = "RdBu")),breaks = seq(-2.5,2.5,length.out = 12),
         show_rownames = F,show_colnames = F)

phm <- bulk$log2[bulk$log2$symbol %in% bulk$deseq2$de90$genesDE,9+which(bulk$info$day == 90)]
pheatmap(mat = phm,scale = "row",color = rev(brewer.pal(n = 11,name = "RdBu")),breaks = seq(-2.5,2.5,length.out = 12),
         show_rownames = F,show_colnames = F,border_color = F)


# Now try with non-LRT version for edgeR

# 150 dpi
rsem <- bulk$rsem[,bulk$info$day == 150]
info <- bulk$info[bulk$info$day == 150,c("genotype","sex")]
rsem_filter <- rsem[rowSums(rsem[,info$genotype == "WT"]) > 0 & rowSums(rsem[,info$genotype == "CKO"]) > 0,]
y <- edgeR::DGEList(counts = rsem_filter,samples = info,group = info$genotype)
y <- edgeR::calcNormFactors(object = y)
design <- model.matrix(object = formula(~ sex + genotype),data = info)
y <- edgeR::estimateDisp(y = y,design = design)
fit <- glmQLFit(y, design)
summary(fit)
de1 <- decideTestsDGE(fit, adjust.method="BH", p.value=0.05)
colnames(design)
