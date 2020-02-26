# Throw everything into UMAP

library(umapr)
library(uwot)

source("./opcBulk/import_opcBulk.R")
source("./opc12/import_opc12.R")
source("./opc90/import_opc90.R")

opc <- cbind(data.frame(row.names = bulk$log2$symbol),
             bulk$log2[,10:ncol(bulk$log2)],
             opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")],
             opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")])
opc <- opc[rowSums(opc) != 0,]


opc.info <- c(paste("bulk",bulk$info$day,bulk$info$genotype,sep = "_"),
              paste("opc12",opc12$info$type[opc12$info$type == "ten-cell"],sep = "_"),
              paste("opc90",opc90$info$type[opc90$info$type == "ten-cell"],sep = "_"))

opc_umap <- uwot::umap(X = t(opc),n_neighbors = 14,learning_rate = 0.5,init = "random")

pdf(file = "./plots/umap_all.pdf",width = 3,height = 3,pointsize = 6)
par(mar = c(4,4,1,1))
plot(x = c(),y = c(),
     xlim = c(-8,16),
     ylim = c(-3,4),
     xlab = "UMAP-1",
     ylab = "UMAP-2",
     xaxs = "i",
     yaxs = "i",
     axes = F)
axis(side = 1,at = seq(-8,16,4))
axis(side = 2,at = seq(-3,4),las = 1)
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
opc_annotation <- data.frame(row.names = names(opc_filter),
                             study = c(rep("bulk",36),rep("opc12",56),rep("opc90",56)),
                             sex = c(as.character(bulk$info$sex),opc12$info$sex[opc12$info$type == "ten-cell"],opc90$info$sex[opc90$info$type == "ten-cell"]),
                             dpi = as.character(c(as.numeric(as.character(bulk$info$day)),rep(12,56),rep(90,56))),
                             genotype = c(as.character(bulk$info$genotype),rep("CKO",112)))

library(pheatmap)
library(RColorBrewer)
png(filename = "./plots/heatmap_all.png",width = 2000,height = 2000,res = 200)
pheatmap(mat = opc_filter,color = rev(brewer.pal(11,"RdBu")[c(1:4,6,8:11)]),scale = "row",show_rownames = F,show_colnames = F,clustering_method = "ward.D2",annotation_col = opc_annotation)
dev.off()
