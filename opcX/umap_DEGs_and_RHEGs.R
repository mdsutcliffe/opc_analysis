# Run UMAP with combined bulk RNA-seq and 10cRNA-seq data

library(uwot)

source("./opc12/import_opc12.R")
source("./opc90/import_opc90.R")
source("./opcBulk/import_opcBulk.R")

source("./opc12/RHEGs_opc12.R")
source("./opc90/RHEGs_opc90.R")

load("./build/bulk_DESeq2_results.RData")

umap_bulk <- cbind(data.frame(row.names = bulk$log2$symbol),
                   bulk$log2[,10:ncol(bulk$log2)])

# Join 12 dpi and 90 dpi 10cRNA-seq datasets
umap_opc12 <- cbind(data.frame(row.names = opc12$log2$symbol[opc12$log2$chr != "ERCC"]),
                    opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")])
umap_opc90 <- cbind(data.frame(row.names = opc12$log2$symbol[opc90$log2$chr != "ERCC"]),
                    opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")])
umap_10c <- cbind(umap_opc12,umap_opc90)

# Z-score the bulk RNA-seq and 10cRNA-seq datasets separately
umap_bulk_scale <- t(scale(t(umap_bulk)))
umap_10c_scale <- t(scale(t(umap_10c)))

# Join scaled bulk RNA-seq and 10cRNA-seq datasets
umap_all_scale <- cbind(umap_bulk_scale,umap_10c_scale)

# List all genes of interest (bulk DEGs + 10c RHEGs)
DEGs_and_RHEGs <- Reduce(union,list(opc12_rheg,
                                    opc90_rheg,
                                    bulk$deseq2$de12$genesDE,
                                    bulk$deseq2$de90$genesDE,
                                    bulk$deseq2$de150$genesDE))

# Retain all genes of interest that are not NA (i.e. SD is zero)
umap_all_scale <- umap_all_scale[DEGs_and_RHEGs,]
umap_all_scale <- umap_all_scale[complete.cases(umap_all_scale),]

# Run UMAP
set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 20)

# Gather information for markers and colors
umap_all.info <- c(paste0("bulk",bulk$info$day,"_",bulk$info$genotype),rep("opc12",sum(opc12$info$type == "ten-cell")),rep("opc90",sum(opc90$info$type == "ten-cell")))
pch_umap <- c(16,1,17,2,15,0,16,16)[as.numeric(factor(umap_all.info))]
col_umap <- c(rep("#000000",6),"#5e3c99","#e66101")[as.numeric(factor(umap_all.info))]

# Plot UMAP of all datasets
pdf(file = "./plots/UMAP_DEGs_and_RHEGs.pdf",width = 2.25,height = 2.25,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     lwd = 0.5,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5,las = 1)
axis(side = 2,lwd = 0.5,las = 1)
legend(x = "topright",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#5e3c99","#e66101"),
       pt.lwd = 0.5,
       box.lwd = 0.5,
       cex = 0.5)
dev.off()

pch_umap <- c(16,1,16,1,16,1,16,16)[as.numeric(factor(umap_all.info))]
col_umap <- c(rep("#000000",6),"#5e3c99","#e66101")[as.numeric(factor(umap_all.info))]

# Plot bulk150
pdf(file = "./plots/UMAP_DEGs_and_RHEGs_bulk150.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
plot(x = embedding_umap[grepl("bulk150",umap_all.info),],
     pch = pch_umap[grepl("bulk150",umap_all.info)],
     col = col_umap[grepl("bulk150",umap_all.info)],
     lwd = 0.5/0.75,
     cex = 0.75,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5/0.75,las = 1)
axis(side = 2,lwd = 0.5/0.75,las = 1)
lines(x = c(median(embedding_umap[umap_all.info == "bulk150_WT",1]),median(embedding_umap[umap_all.info == "bulk150_CKO",1])),
      y = c(median(embedding_umap[umap_all.info == "bulk150_WT",2]),median(embedding_umap[umap_all.info == "bulk150_CKO",2])),
      lwd = 0.5/0.75)

# legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),col = "#000000",pt.lwd = 0.5,box.lwd = 0.5,cex = 6/7)
text(x = c(-2,-0.5),y = c(-1,2.5),labels = c("Control","N1P"),cex = 6/7)
dev.off()

# Plot bulk12
pdf(file = "./plots/UMAP_DEGs_and_RHEGs_bulk12.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
plot(x = embedding_umap[grepl("bulk12",umap_all.info),],
     pch = pch_umap[grepl("bulk12",umap_all.info)],
     col = col_umap[grepl("bulk12",umap_all.info)],
     lwd = 0.5/0.75,
     cex = 0.75,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5/0.75,las = 1)
axis(side = 2,lwd = 0.5/0.75,las = 1)
lines(x = c(median(embedding_umap[umap_all.info == "bulk12_WT",1]),median(embedding_umap[umap_all.info == "bulk12_CKO",1])),
      y = c(median(embedding_umap[umap_all.info == "bulk12_WT",2]),median(embedding_umap[umap_all.info == "bulk12_CKO",2])),
      lwd = 0.5/0.75)
# legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),col = "#000000",pt.lwd = 0.5,box.lwd = 0.5,cex = 6/7)
text(x = c(-2,2),y = c(0.5,-1),labels = c("Control","N1P"),cex = 6/7)
dev.off()

# Perform k-means clustering on bulk90 N1P samples
cluster90 <- kmeans(x = embedding_umap[umap_all.info == "bulk90_CKO",],centers = 2)

# Plot bulk90
pdf(file = "./plots/UMAP_DEGs_and_RHEGs_bulk90.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
plot(x = embedding_umap[grepl("bulk90",umap_all.info),],
     pch = pch_umap[grepl("bulk90",umap_all.info)],
     col = col_umap[grepl("bulk90",umap_all.info)],
     lwd = 0.5/0.75,
     cex = 0.75,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5/0.75,las = 1)
axis(side = 2,lwd = 0.5/0.75,las = 1)
lines(x = c(median(embedding_umap[umap_all.info == "bulk90_WT",1]),median(embedding_umap[umap_all.info == "bulk90_CKO",1][cluster90$cluster == 1])),
      y = c(median(embedding_umap[umap_all.info == "bulk90_WT",2]),median(embedding_umap[umap_all.info == "bulk90_CKO",2][cluster90$cluster == 1])),
      lwd = 0.5/0.75)
lines(x = c(median(embedding_umap[umap_all.info == "bulk90_WT",1]),median(embedding_umap[umap_all.info == "bulk90_CKO",1][cluster90$cluster == 2])),
      y = c(median(embedding_umap[umap_all.info == "bulk90_WT",2]),median(embedding_umap[umap_all.info == "bulk90_CKO",2][cluster90$cluster == 2])),
      lwd = 0.5/0.75)
# legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),col = "#000000",pt.lwd = 0.5,box.lwd = 0.5,cex = 6/7)
text(x = c(-1,2,-1),y = c(-1.5,-0.25,1),labels = c("Control","N1P","N1P"),cex = 6/7)
dev.off()

library(e1071)
dat <- data.frame(embedding_umap[grepl("opc",umap_all.info),],y = factor(c(rep("12 dpi",56),rep("90 dpi",56))))
svmfit <- svm(y ~ .,data = dat,kernel = "linear",scale = F)
beta = drop(t(svmfit$coefs)%*%embedding_umap[grepl("opc",umap_all.info),][svmfit$index,])
beta0 = svmfit$rho

# Plot 10cRNA-seq
pdf(file = "./plots/UMAP_DEGs_and_RHEGs_10c.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
plot(x = embedding_umap[grepl("opc",umap_all.info),],
     pch = pch_umap[grepl("opc",umap_all.info)],
     col = col_umap[grepl("opc",umap_all.info)],
     lwd = 0.5/0.75,
     cex = 0.75,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
abline(beta0 / beta[2], -beta[1] / beta[2],lwd = 0.5/0.75)
axis(side = 1,lwd = 0.5/0.75,las = 1)
axis(side = 2,lwd = 0.5/0.75,las = 1)
# legend(x = "topright",legend = c("10c-12","10c-90"),pch = c(1,16),col = "#000000",pt.lwd = 0.5,box.lwd = 0.5,cex = 6/7)
text(x = c(-1.5,2),y = c(-1.5,1),labels = c("12 dpi","90 dpi"),cex = 6/7)
dev.off()



pdf(file = "./plots/UMAP_DEGs_and_RHEGs_12_panelE.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
plot(x = embedding_umap[grepl("bulk12",umap_all.info) | grepl("opc12",umap_all.info),],
     pch = pch_umap[grepl("bulk12",umap_all.info) | grepl("opc12",umap_all.info)],
     col = col_umap[grepl("bulk12",umap_all.info) | grepl("opc12",umap_all.info)],
     lwd = 0.5/0.75,
     cex = 0.75,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5/0.75,las = 1)
axis(side = 2,lwd = 0.5/0.75,las = 1)
lines(x = c(median(embedding_umap[umap_all.info == "bulk12_WT",1]),median(embedding_umap[umap_all.info == "bulk12_CKO",1])),
      y = c(median(embedding_umap[umap_all.info == "bulk12_WT",2]),median(embedding_umap[umap_all.info == "bulk12_CKO",2])),
      lwd = 0.5/0.75)
# text(x = c(-2,2),y = c(0.5,-1),labels = c("Control","N1P"),cex = 6/7)
dev.off()


pdf(file = "./plots/UMAP_DEGs_and_RHEGs_90_panelF.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F,family = "ArialMT")
par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
plot(x = embedding_umap[grepl("bulk90",umap_all.info) | grepl("opc90",umap_all.info),],
     pch = pch_umap[grepl("bulk90",umap_all.info) | grepl("opc90",umap_all.info)],
     col = col_umap[grepl("bulk90",umap_all.info) | grepl("opc90",umap_all.info)],
     lwd = 0.5/0.75,
     cex = 0.75,
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     frame = F,
     axes = F,
     xaxs = "i",
     yaxs = "i")
axis(side = 1,lwd = 0.5/0.75,las = 1)
axis(side = 2,lwd = 0.5/0.75,las = 1)
lines(x = c(median(embedding_umap[umap_all.info == "bulk90_WT",1]),median(embedding_umap[umap_all.info == "bulk90_CKO",1][cluster90$cluster == 1])),
      y = c(median(embedding_umap[umap_all.info == "bulk90_WT",2]),median(embedding_umap[umap_all.info == "bulk90_CKO",2][cluster90$cluster == 1])),
      lwd = 0.5/0.75)
lines(x = c(median(embedding_umap[umap_all.info == "bulk90_WT",1]),median(embedding_umap[umap_all.info == "bulk90_CKO",1][cluster90$cluster == 2])),
      y = c(median(embedding_umap[umap_all.info == "bulk90_WT",2]),median(embedding_umap[umap_all.info == "bulk90_CKO",2][cluster90$cluster == 2])),
      lwd = 0.5/0.75)
# text(x = c(-2,2),y = c(0.5,-1),labels = c("Control","N1P"),cex = 6/7)
dev.off()