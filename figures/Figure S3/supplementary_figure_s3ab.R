# Supplementary Figure S3 - panels A and B

library(pheatmap)
library(RColorBrewer)

source("./functions/signatureMatrix.R")
source("./import/import_opcBulk.R")

signature <- signatureMatrix(geneList = bulk$tpm$symbol)

sig_mat <- t(as.matrix(signature$matrix[,2:ncol(signature$matrix)]))
rownames(sig_mat)[rownames(sig_mat) %in% "MO"] <- "Myelinating\noligodendrocyte"

pheatmap(mat = sig_mat,
         color = rev(brewer.pal(11,"RdBu")),
         scale = "column",
         cellwidth = 0.5,
         cellheight = 20,
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F,
         fontsize = 7,
         annotation = data.frame(row.names = colnames(sig_mat),
                                 group = rep(x = 1:nrow(sig_mat),each = 40)),
         annotation_colors = list(group = rev(brewer.pal(6,"Set1"))),
         annotation_legend = F,
         annotation_names_col = F,
         filename = "./Figure S3/figure_s3a.pdf",
         width = 3.125,
         height = 2.34375)

write.table(x = signature$tpm,file = "./Figure S3/mixture_signature.txt",quote = F,sep = "\t",row.names = F)
write.table(x = signature$matrix,file = "./Figure S3/signature.txt",quote = F,sep = "\t",row.names = F)

# Run cibersort, get results, and place in external folder

f.cibersort <- "./Figure S3/CIBERSORTx_signature_matrix_validation.txt"

cibersort_signature <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort_signature) <- cibersort_signature$Mixture
row.names(cibersort_signature)[row.names(cibersort_signature) == "MO"] <- "Myelinating\noligodendrocyte"
names(cibersort_signature)[names(cibersort_signature) == "MO"] <- "Myelinating\noligodendrocyte"
cibersort_signature <- cibersort_signature[,2:(which(names(cibersort_signature) == "P.value") - 1)] / cibersort_signature$Absolute.score..sig.score.

pdf(file = "./Figure S3/figure_s3b.pdf",width = 3.125,height = 2.34375,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.75,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
barplot(t(as.matrix(cibersort_signature[nrow(cibersort_signature):1,ncol(cibersort_signature):1]))*100,
        horiz = T,
        las = 1,
        lwd = 0.5/0.75,
        col = brewer.pal(6,"Set1"),
        xlab = "Relative proportion")
dev.off()

source("./functions/signature_matrix_pericytes.R")
source("./import/import_opcBulk.R")

signature <- signatureMatrix(geneList = bulk$tpm$symbol)