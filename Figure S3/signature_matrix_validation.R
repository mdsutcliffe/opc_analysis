# Signature matrix validation by cibersort

source("./functions/signatureMatrix.R")
source("./opcBulk/import_opcBulk.R")

signature <- signatureMatrix(geneList = bulk$tpm$symbol)

sig_mat <- t(as.matrix(signature$matrix[,2:ncol(signature$matrix)]))
rownames(sig_mat)[rownames(sig_mat) %in% "MO"] <- "Myelinating\noligodendrocyte"

pdf(file = "./plots/figure_s3a.pdf",width = 3.5,height = 2.25,family = "ArialMT",pointsize = 7)
pheatmap(mat = sig_mat,
         color = rev(brewer.pal(11,"RdBu")),
         scale = "column",
         cellwidth = 0.5,
         cellheight = 20,
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F,
         fontsize = 7)
dev.off()

write.table(x = signature$tpm,file = "./Figure S3/mixture_signature.txt",quote = F,sep = "\t",row.names = F)
write.table(x = signature$matrix,file = "./Figure S3/signature.txt",quote = F,sep = "\t",row.names = F)

# Run cibersort, get results, and place in external folder

f.cibersort <- "./Figure S3/CIBERSORTx_signature_matrix_validation.txt"

cibersort_signature <- read.table(file = f.cibersort,header = T,sep = "\t")
row.names(cibersort_signature) <- cibersort_signature$Mixture
row.names(cibersort_signature)[row.names(cibersort_signature) == "MO"] <- "Myelinating\noligodendrocyte"
names(cibersort_signature)[names(cibersort_signature) == "MO"] <- "Myelinating\noligodendrocyte"
cibersort_signature <- cibersort_signature[,2:(which(names(cibersort_signature) == "P.value") - 1)] / cibersort_signature$Absolute.score..sig.score.

pdf(file = "./plots/figure_s3b.pdf",width = 2.25,height = 2,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
barplot(t(as.matrix(cibersort_signature[nrow(cibersort_signature):1,ncol(cibersort_signature):1]))*100,
        horiz = T,
        yaxs = "i",
        las = 1,
        lwd = 0.5/0.75,
        col = brewer.pal(6,"Set1"),
        xlab = "Relative proportion")
axis(side = 2,at = seq(0,8,1.25),labels = F,lwd = 0.5/0.75)
dev.off()

source("./functions/signature_matrix_pericytes.R")
source("./opcBulk/import_opcBulk.R")

signature <- signatureMatrix(geneList = bulk$tpm$symbol)