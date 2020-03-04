# Throw everything into UMAP

# library(umapr)
library(uwot)

source("./opcBulk/import_opcBulk.R")
source("./opc12/import_opc12.R")
source("./opc90/import_opc90.R")
source("./opcBulk/cleancounts_bulk.R")

opc12F <- read.table(file = "./data/opc12_true/cleancounts_female_true.tsv")
opc12M <- read.table(file = "./data/opc12_true/cleancounts_male_true.tsv")
opc90F <- read.table(file = "./data/opc90_true/cleancounts_female_true.tsv")
opc90M <- read.table(file = "./data/opc90_true/cleancounts_male_true.tsv")

allCommonGenes <- Reduce(intersect,list(row.names(opc12F),row.names(opc12M),
                      row.names(opc90F),row.names(opc90M),
                      row.names(bulk_clean$bulk12),row.names(bulk_clean$bulk90),row.names(bulk_clean$bulk150)))

opc <- cbind(data.frame(row.names = bulk$log2$symbol),
             bulk$log2[,10:ncol(bulk$log2)],
             opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")],
             opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")])
opc <- opc[allCommonGenes,]

opc.info <- c(paste("bulk",bulk$info$day,bulk$info$genotype,sep = "_"),
              paste("opc12",opc12$info$type[opc12$info$type == "ten-cell"],sep = "_"),
              paste("opc90",opc90$info$type[opc90$info$type == "ten-cell"],sep = "_"))

opc_umap <- uwot::umap(X = t(opc),n_neighbors = 14)
plot(opc_umap)

pdf(file = "./plots/umap_all.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-6,6),
     ylim = c(-4,8),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1,at = seq(-6,6,2))
axis(side = 2,at = seq(-4,8,2),las = 1)
points(opc_umap[grepl(pattern = "WT",x = opc.info),],pch = 1,lwd = 0.5)
points(opc_umap[grepl(pattern = "CKO",x = opc.info),],pch = 16)

points(opc_umap[grepl(pattern = "opc12",x = opc.info),],pch = 16,col = "#e41a1c")
points(opc_umap[grepl(pattern = "opc90",x = opc.info),],pch = 16,col = "#377eb8")
legend(x = "topright",legend = c("bulk WT (all days)","bulk CKO (all days)","10c-12dpi","10c-90dpi"),col = c("#000000","#000000","#e41a1c","#377eb8"),pch = c(1,16,16,16),pt.lwd = c(0.5,0.5,0.5,0.5),box.lwd = 0.5)
dev.off()

opc_pca <- prcomp(x = t(opc))
library(ggbiplot)
pca_plot <- ggbiplot(opc_pca,var.axes = F,groups = factor(opc.info)) + theme_bw()
ggsave(filename = "./plots/pca_all.pdf",plot = pca_plot,width = 6,height = 6,units = "in")

opc_filter <- opc[rowSums(opc > 0) > 37,]
opc_annotation <- data.frame(row.names = names(opc),
                             study = c(rep("bulk",36),rep("opc12",56),rep("opc90",56)),
                             sex = c(as.character(bulk$info$sex),opc12$info$sex[opc12$info$type == "ten-cell"],opc90$info$sex[opc90$info$type == "ten-cell"]),
                             dpi = as.character(c(as.numeric(as.character(bulk$info$day)),rep(12,56),rep(90,56))),
                             genotype = c(as.character(bulk$info$genotype),rep("CKO",112)))

library(pheatmap)
library(RColorBrewer)
png(filename = "./plots/heatmap_all.png",width = 2000,height = 2000,res = 200)
pheatmap(mat = opc,color = rev(brewer.pal(11,"RdBu")[c(1:4,6,8:11)]),scale = "row",show_rownames = F,show_colnames = F,clustering_method = "ward.D2",annotation_col = opc_annotation)
dev.off()



bulk_only <- bulk$log2[bulk$log$symbol %in% Reduce(intersect,list(row.names(bulk_clean$bulk12),row.names(bulk_clean$bulk90),row.names(bulk_clean$bulk150))),10:ncol(bulk$log2)]
bulk_umap <- uwot::umap(X = t(bulk_only),n_neighbors = 14)
plot(bulk_umap)
pdf(file = "./plots/umap_bulk.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-4,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(bulk_umap[bulk$info$day == 12 & bulk$info$genotype == "WT",],pch = 1,lwd = 0.5,col = "#3182bd")
points(bulk_umap[bulk$info$day == 12 & bulk$info$genotype == "CKO",],pch = 16,lwd = 0.5,col = "#3182bd")
points(bulk_umap[bulk$info$day == 90 & bulk$info$genotype == "WT",],pch = 1,lwd = 0.5,col = "#31a354")
points(bulk_umap[bulk$info$day == 90 & bulk$info$genotype == "CKO",],pch = 16,lwd = 0.5,col = "#31a354")
points(bulk_umap[bulk$info$day == 150 & bulk$info$genotype == "WT",],pch = 1,lwd = 0.5,col = "#de2d26")
points(bulk_umap[bulk$info$day == 150 & bulk$info$genotype == "CKO",],pch = 16,lwd = 0.5,col = "#de2d26")

legend(x = "topright",
       legend = c("12 dpi WT","12 dpi N1P","90 dpi WT","90 dpi N1P","150 dpi WT","150 dpi N1P"),
       col = rep(x = c("#3182bd","#31a354","#de2d26"),each = 2),
       pch = c(1,16),pt.lwd = c(0.5,0.5,0.5,0.5),box.lwd = 0.5)
dev.off()


tencell <- cbind(opc12$log2[,9+which(opc12$info$type == "ten-cell")],opc90$log2[,9+which(opc90$info$type == "ten-cell")])
tencell <- tencell[opc12$log2$symbol %in% Reduce(intersect,list(row.names(opc12F),row.names(opc12M),row.names(opc90F),row.names(opc90M))),]
tencell_umap <- uwot::umap(X = t(tencell),n_neighbors = 14)
plot(tencell_umap)

pdf(file = "./plots/umap_10c.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-4,6),
     ylim = c(-2,2),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(tencell_umap[(1:56)[opc12$info$sex == "female"],],pch = 1,lwd = 2,col = "#3182bd")
points(tencell_umap[(1:56)[opc12$info$sex == "male"],],pch = 4,lwd = 2,col = "#3182bd")
points(tencell_umap[(57:112)[opc90$info$sex == "female"],],pch = 1,lwd = 2,col = "#31a354")
points(tencell_umap[(57:112)[opc90$info$sex == "male"],],pch = 4,lwd = 2,col = "#31a354")
legend(x = "topright",
       legend = c("12 dpi Female","12 dpi Male","90 dpi Female","90 dpi Male"),
       col = rep(x = c("#3182bd","#31a354"),each = 2),
       pch = rep(x = c(1,4),2),
       pt.lwd = 2, box.lwd = 0.5)
dev.off()


opc12_10c <- opc12$log2[opc12$log2$symbol %in% intersect(row.names(opc12F),row.names(opc12M)),c(3,9+which(opc12$info$type == "ten-cell"))]
# opc12_10c <- opc12_10c[opc12$log,]
set.seed(12)

tempumap <- t(opc12_10c)
colnames(tempumap) <- tempumap[1,]
tempumap <- apply(tempumap,2,as.numeric)
tempumap <- as.matrix(tempumap[2:nrow(tempumap),])

opc12_10c_umap <- uwot::umap(X = tempumap,n_neighbors = 7,)
opc12_10c_umap <- umapr::umap(data = tempumap,n_neighbors = 7)
plot(opc12_10c_umap)
umapr::run_umap_shiny(cbind(tempumap[,c("Pdgfra")],opc12_10c_umap))

xx <- cbind(tempumap[,"Pdgfra"],opc12_10c_umap)
umapr::run_umap_shiny(xx)

yy <- umapr::umap(data = tempumap,n_neighbors = 7)
umapr::run_umap_shiny(yy)
embedding <- yy
umap1_cors <- apply(embedding[,1:(ncol(embedding)-2)],2,function(x) cor(x,embedding[,(ncol(embedding)-1)]))
umap2_cors <- apply(embedding[,1:(ncol(embedding)-2)],2,function(x) cor(x,embedding[,(ncol(embedding))]))

pdf(file = "./plots/umap_10c_12dpi.pdf",width = 2.25,height = 2,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1),mgp = c(2,1,0))
plot(x = c(),y = c(),
     xlim = c(-6,6),
     ylim = c(-8,8),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1,at = seq(-8,8,4))
points(opc12_10c_umap[opc12$info$sex[opc12$info$type == "ten-cell"] == "female",],pch = 1,lwd = 2,col = "#3182bd")
points(opc12_10c_umap[opc12$info$sex[opc12$info$type == "ten-cell"] == "male",],pch = 4,lwd = 2,col = "#3182bd")
legend(x = "topright",
       legend = c("12 dpi Female","12 dpi Male"),
       col = rep(x = "#3182bd",2),
       pch = c(1,4),
       pt.lwd = 2, box.lwd = 0.5)
dev.off()

library(RColorBrewer)
cb <- brewer.pal(11,"RdBu")
colPal <- colorRampPalette(colors = c())
pdf(file = "./plots/umap_10c_12dpi.pdf",width = 2.25,height = 2,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1),mgp = c(2,1,0))
plot(x = c(),y = c(),
     xlim = c(-12,6),
     ylim = c(-8,8),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1,at = seq(-8,8,4))
points(opc12_10c_umap[opc12$info$sex[opc12$info$type == "ten-cell"] == "female",],pch = 1,lwd = 2,col = "#3182bd")
points(opc12_10c_umap[opc12$info$sex[opc12$info$type == "ten-cell"] == "male",],pch = 4,lwd = 2,col = "#3182bd")

points(opc12_10c_umap,
       pch = c(1,4)[as.numeric(opc12$info$sex[opc12$info$type == "ten-cell"] == "female") + 1],
       col = cb[as.numeric(cut(as.numeric(opc12_10c[opc12_10c$symbol == "Epha5",2:ncol(opc12_10c)]),breaks = 11))],
       lwd = 4)
legend(x = "topleft",legend = c("Low","","High"),col = cb[c(1,6,11)],fill = T)


legend(x = "topright",
       legend = c("12 dpi Female","12 dpi Male"),
       col = rep(x = "#3182bd",2),
       pch = c(1,4),
       pt.lwd = 2, box.lwd = 0.5)
dev.off()


View(umap1_cors * umap2_cors)

opc90_10c <- opc90$log2[opc90$log2$symbol %in% intersect(row.names(opc90F),row.names(opc90M)),9+which(opc90$info$type == "ten-cell")]
# opc90_10c <- opc90_10c[rowSums(opc90_10c > 0) > 7,]
opc90_10c_umap <- uwot::umap(X = t(opc90_10c),n_neighbors = 7)
plot(opc90_10c_umap)
pdf(file = "./plots/umap_10c_90dpi.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-4,4),
     ylim = c(-2,2),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(opc90_10c_umap[opc90$info$sex[opc90$info$type == "ten-cell"] == "female",],pch = 1,lwd = 2,col = "#31a354")
points(opc90_10c_umap[opc90$info$sex[opc90$info$type == "ten-cell"] == "male",],pch = 4,lwd = 2,col = "#31a354")
legend(x = "topright",
       legend = c("90 dpi Female","90 dpi Male"),
       col = rep(x = "#31a354",2),
       pch = c(1,4),
       pt.lwd = 2, box.lwd = 0.5)
dev.off()

uwot::



# Bulk using only DE genes -----
source("./opcBulk/DESeq2_opcBulk.R")
bulk_only <- bulk$log2[bulk$log$symbol %in% Reduce(f = union,x = list(res_plus_12$genesDE,
                                                                      res_plus_90$genesDE,
                                                                      res_plus_150$genesDE)),10:ncol(bulk$log2)]
set.seed(0)
bulk_umap <- uwot::umap(X = t(bulk_only),n_neighbors = 14)
plot(bulk_umap)
pdf(file = "./plots/umap_bulk_DEgenes.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-10,4),
     ylim = c(-4,6),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(bulk_umap[bulk$info$day == 12 & bulk$info$genotype == "WT",],pch = 1,lwd = 0.5,col = "#3182bd")
points(bulk_umap[bulk$info$day == 12 & bulk$info$genotype == "CKO",],pch = 16,lwd = 0.5,col = "#3182bd")
points(bulk_umap[bulk$info$day == 90 & bulk$info$genotype == "WT",],pch = 1,lwd = 0.5,col = "#31a354")
points(bulk_umap[bulk$info$day == 90 & bulk$info$genotype == "CKO",],pch = 16,lwd = 0.5,col = "#31a354")
points(bulk_umap[bulk$info$day == 150 & bulk$info$genotype == "WT",],pch = 1,lwd = 0.5,col = "#de2d26")
points(bulk_umap[bulk$info$day == 150 & bulk$info$genotype == "CKO",],pch = 16,lwd = 0.5,col = "#de2d26")
legend(x = "topright",
       legend = c("12 dpi WT","12 dpi N1P","90 dpi WT","90 dpi N1P","150 dpi WT","150 dpi N1P"),
       col = rep(x = c("#3182bd","#31a354","#de2d26"),each = 2),
       pch = c(1,16),pt.lwd = c(0.5,0.5,0.5,0.5),box.lwd = 0.5)
dev.off()

# Ten-cell using only the RHEGs -----
source("./opc12/RHEGs_opc12.R")
source("./opc90/RHEGs_opc90.R")
tencell <- cbind(data.frame(row.names = opc12$log2$symbol),
                 opc12$log2[,9+which(opc12$info$type == "ten-cell")],
                 opc90$log2[,9+which(opc90$info$type == "ten-cell")])
tencell <- tencell[Reduce(f = union,x = list(opc12_rheg,opc90_rheg)),]


set.seed(0)
tencell_umap <- uwot::umap(X = t(tencell),n_neighbors = 7)
plot(tencell_umap)



pdf(file = "./plots/umap_10c.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-4,4),
     ylim = c(-4,4),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(tencell_umap[(1:56)[opc12$info$sex == "female"],],pch = 1,lwd = 2,col = "#3182bd")
points(tencell_umap[(1:56)[opc12$info$sex == "male"],],pch = 4,lwd = 2,col = "#3182bd")
points(tencell_umap[(57:112)[opc90$info$sex == "female"],],pch = 1,lwd = 2,col = "#31a354")
points(tencell_umap[(57:112)[opc90$info$sex == "male"],],pch = 4,lwd = 2,col = "#31a354")
legend(x = "topright",
       legend = c("12 dpi Female","12 dpi Male","90 dpi Female","90 dpi Male"),
       col = rep(x = c("#3182bd","#31a354"),each = 2),
       pch = rep(x = c(1,4),2),
       pt.lwd = 2, box.lwd = 0.5)
dev.off()


# opc12 RHEGs only -----
opc12_10c <- opc12$log2[opc12$log2$symbol %in% opc12_rheg,9+which(opc12$info$type == "ten-cell")]
# opc12_10c <- opc12$log2[opc12$log2$symbol %in% Reduce(f = union,x = list(opc12_rheg,opc12_uniqueF,opc12_uniqueM)),9+which(opc12$info$type == "ten-cell")]
# opc12_10c <- opc12_10c[opc12$log,]
set.seed(0)
opc12_10c_umap <- uwot::umap(X = t(opc12_10c),n_neighbors = 7)
plot(opc12_10c_umap)
pdf(file = "./plots/umap_10c_12dpi.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(opc12_10c_umap[opc12$info$sex[opc12$info$type == "ten-cell"] == "female",],pch = 1,lwd = 2,col = "#3182bd")
points(opc12_10c_umap[opc12$info$sex[opc12$info$type == "ten-cell"] == "male",],pch = 4,lwd = 2,col = "#3182bd")
legend(x = "topright",
       legend = c("12 dpi Female","12 dpi Male"),
       col = rep(x = "#3182bd",2),
       pch = c(1,4),
       pt.lwd = 2, box.lwd = 0.5)
dev.off()

# opc90 RHEGs only -----
opc90_10c <- opc90$log2[opc90$log2$symbol %in% opc90_rheg,9+which(opc90$info$type == "ten-cell")]
# opc90_10c <- opc90$log2[opc90$log2$symbol %in% Reduce(f = union,x = list(opc90_rheg,opc90_uniqueF,opc90_uniqueM)),9+which(opc90$info$type == "ten-cell")]
# opc90_10c <- opc90_10c[rowSums(opc90_10c > 0) > 7,]
set.seed(0)
opc90_10c_umap <- uwot::umap(X = t(opc90_10c),n_neighbors = 7)
plot(opc90_10c_umap)
pdf(file = "./plots/umap_10c_90dpi.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-3,3),
     ylim = c(-3,3),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(opc90_10c_umap[opc90$info$sex[opc90$info$type == "ten-cell"] == "female",],pch = 1,lwd = 2,col = "#31a354")
points(opc90_10c_umap[opc90$info$sex[opc90$info$type == "ten-cell"] == "male",],pch = 4,lwd = 2,col = "#31a354")
legend(x = "topright",
       legend = c("90 dpi Female","90 dpi Male"),
       col = rep(x = "#31a354",2),
       pch = c(1,4),
       pt.lwd = 2, box.lwd = 0.5)
dev.off()

# Everything together, only DE + RHEGs

allCommonGenes <- Reduce(union,list(opc12_rheg,
                                        opc90_rheg,
                                        res_plus_12$genesDE,
                                        res_plus_90$genesDE,
                                        res_plus_150$genesDE))

opc <- cbind(data.frame(row.names = bulk$log2$symbol),
             bulk$log2[,10:ncol(bulk$log2)],
             opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")],
             opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")])
opc <- opc[allCommonGenes,]

opc.info <- c(paste("bulk",bulk$info$day,bulk$info$genotype,sep = "_"),
              paste("opc12",opc12$info$type[opc12$info$type == "ten-cell"],sep = "_"),
              paste("opc90",opc90$info$type[opc90$info$type == "ten-cell"],sep = "_"))

set.seed(0)
opc_umap <- uwot::umap(X = t(opc),n_neighbors = 14)
plot(opc_umap)

pdf(file = "./plots/umap_all_de_rheg.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-4,10),
     ylim = c(-10,6),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(opc_umap[grepl(pattern = "WT",x = opc.info),],pch = 1,lwd = 0.5)
points(opc_umap[grepl(pattern = "CKO",x = opc.info),],pch = 16)

points(opc_umap[grepl(pattern = "opc12",x = opc.info),],pch = 16,col = "#e41a1c")
points(opc_umap[grepl(pattern = "opc90",x = opc.info),],pch = 16,col = "#377eb8")
legend(x = "topright",legend = c("bulk WT (all days)","bulk CKO (all days)","10c-12dpi","10c-90dpi"),col = c("#000000","#000000","#e41a1c","#377eb8"),pch = c(1,16,16,16),pt.lwd = c(0.5,0.5,0.5,0.5),box.lwd = 0.5)
dev.off()

# Everything together, DE + union of M/F unique/rhegs
source("./opcBulk/runDESeq2_opcBulk.R")
allCommonGenes <- Reduce(union,list(opc12_rheg,opc12_uniqueF,opc12_uniqueM,
                                    opc90_rheg,opc90_uniqueF,opc90_uniqueM,
                                    bulk$deseq2$de12$genesDE,
                                    bulk$deseq2$de90$genesDE,
                                    bulk$deseq2$de150$genesDE))

scale_bulk <- t(scale(t(cbind(data.frame(row.names = bulk$log2$symbol),bulk$log2[,10:ncol(bulk$log2)]))))
scale_10c <- t(scale(t(cbind(data.frame(row.names = bulk$log2$symbol),opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")],opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")]))))

opc_scale <- cbind(scale_bulk,scale_10c)
opc_scale <- opc_scale[allCommonGenes,]
opc_scale <- opc_scale[complete.cases(opc_scale),]


opc <- cbind(data.frame(row.names = bulk$log2$symbol),
             bulk$log2[,10:ncol(bulk$log2)],
             opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")],
             opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")])
opc <- opc[allCommonGenes,]

opc.info <- c(paste("bulk",bulk$info$day,bulk$info$genotype,sep = "_"),
              paste("opc12",opc12$info$type[opc12$info$type == "ten-cell"],sep = "_"),
              paste("opc90",opc90$info$type[opc90$info$type == "ten-cell"],sep = "_"))

set.seed(0)
opc_umap <- uwot::umap(X = t(opc_scale),n_neighbors = 12)
plot(opc_umap)

pdf(file = "./plots/umap_all_de_rheg_unionMF_scale.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-2,2),
     ylim = c(-3,4),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(opc_umap[grepl(pattern = "WT",x = opc.info),],pch = 1,lwd = 0.5)
points(opc_umap[grepl(pattern = "CKO",x = opc.info),],pch = 16)
points(opc_umap[grepl(pattern = "150_CKO",x = opc.info),],pch = 17)

points(opc_umap[grepl(pattern = "opc12",x = opc.info),],pch = 16,col = "#e41a1c")
points(opc_umap[grepl(pattern = "opc90",x = opc.info),],pch = 16,col = "#377eb8")
legend(x = "topright",legend = c("bulk WT (all days)","bulk CKO (all days)","10c-12dpi","10c-90dpi"),col = c("#000000","#000000","#e41a1c","#377eb8"),pch = c(1,16,16,16),pt.lwd = c(0.5,0.5,0.5,0.5),box.lwd = 0.5)
dev.off()






allCommonGenes <- Reduce(union,list(opc12_rheg,
                                    opc90_rheg,
                                    bulk$deseq2$de12$genesDE,
                                    bulk$deseq2$de90$genesDE,
                                    bulk$deseq2$de150$genesDE))

scale_bulk <- t(scale(t(cbind(data.frame(row.names = bulk$log2$symbol),bulk$log2[,10:ncol(bulk$log2)]))))
scale_10c <- t(scale(t(cbind(data.frame(row.names = bulk$log2$symbol),opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")],opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")]))))

opc_scale <- cbind(scale_bulk,scale_10c)
opc_scale <- opc_scale[allCommonGenes,]
opc_scale <- opc_scale[complete.cases(opc_scale),]


opc <- cbind(data.frame(row.names = bulk$log2$symbol),
             bulk$log2[,10:ncol(bulk$log2)],
             opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")],
             opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")])
opc <- opc[allCommonGenes,]

opc.info <- c(paste("bulk",bulk$info$day,bulk$info$genotype,sep = "_"),
              paste("opc12",opc12$info$type[opc12$info$type == "ten-cell"],sep = "_"),
              paste("opc90",opc90$info$type[opc90$info$type == "ten-cell"],sep = "_"))

set.seed(0)
opc_umap <- uwot::umap(X = t(opc_scale),n_neighbors = 12)
plot(opc_umap)

pdf(file = "./plots/umap_all_de_rhegsONLY_unionMF_scale.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-2,2),
     ylim = c(-3,4),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1)
axis(side = 2,las = 1)
points(opc_umap[grepl(pattern = "WT",x = opc.info),],pch = 1,lwd = 0.5)
points(opc_umap[grepl(pattern = "CKO",x = opc.info),],pch = 16)
points(opc_umap[grepl(pattern = "150_CKO",x = opc.info),],pch = 17)

points(opc_umap[grepl(pattern = "opc12",x = opc.info),],pch = 16,col = "#e41a1c")
points(opc_umap[grepl(pattern = "opc90",x = opc.info),],pch = 16,col = "#377eb8")
legend(x = "topright",legend = c("bulk WT (all days)","bulk CKO (all days)","10c-12dpi","10c-90dpi"),col = c("#000000","#000000","#e41a1c","#377eb8"),pch = c(1,16,16,16),pt.lwd = c(0.5,0.5,0.5,0.5),box.lwd = 0.5)
dev.off()

pch_scale <- sapply(opc.info,function(x) switch(x,"opc12_ten-cell" = 16,"opc90_ten-cell" = 16,"bulk_12_WT" = 1,"bulk_90_WT" = 1,"bulk_150_WT" = 1,"bulk_12_CKO" = 16,"bulk_90_CKO" = 18,"bulk_150_CKO" = 8))
col_scale <- sapply(opc.info,function(x) switch(x,"opc12_ten-cell" = "#e41a1c","opc90_ten-cell" = "#377eb8","bulk_12_WT" = "#000000","bulk_90_WT" = "#000000","bulk_150_WT" = "#000000","bulk_12_CKO" = "#000000","bulk_90_CKO" = "#000000","bulk_150_CKO" = "#000000"))
pdf(file = "./plots/umap_all_de_rhegsONLY_unionMF_scale.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = opc_umap,
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     # xaxs = "i",
     # yaxs = "i",
     axes = F,
     pch = pch_scale,
     col = col_scale)
axis(side = 1)
axis(side = 2,las = 1)
legend(x = "topright",legend = c("bulk WT (all days)","bulk CKO (all days)","10c-12dpi","10c-90dpi"),col = c("#000000","#000000","#e41a1c","#377eb8"),pch = c(1,16,16,16),pt.lwd = c(0.5,0.5,0.5,0.5),box.lwd = 0.5)
dev.off()