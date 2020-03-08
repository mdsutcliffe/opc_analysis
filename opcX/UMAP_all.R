setwd("~")
library(uwot)

setwd("~/Github/opc_analysis/")
source("./opc12/import_opc12.R")
source("./opc90/import_opc90.R")
source("./opcBulk/import_opcBulk.R")

source("./opc12/RHEGs_opc12.R")
source("./opc90/RHEGs_opc90.R")
source("./opcBulk/runDESeq2_opcBulk.R")

umap_bulk <- cbind(data.frame(row.names = bulk$log2$symbol),
                   bulk$log2[,10:ncol(bulk$log2)])
umap_opc12 <- cbind(data.frame(row.names = opc12$log2$symbol[opc12$log2$chr != "ERCC"]),
                    opc12$log2[opc12$log2$chr != "ERCC",9 + which(opc12$info$type == "ten-cell")])
umap_opc90 <- cbind(data.frame(row.names = opc12$log2$symbol[opc90$log2$chr != "ERCC"]),
                    opc90$log2[opc90$log2$chr != "ERCC",9 + which(opc90$info$type == "ten-cell")])
umap_10c <- cbind(umap_opc12,umap_opc90)

umap_bulk_scale <- t(scale(t(umap_bulk)))
umap_10c_scale <- t(scale(t(umap_10c)))
umap_all_scale <- cbind(umap_bulk_scale,umap_10c_scale)

genes_DE_and_het <- Reduce(union,list(opc12_rheg,opc12_uniqueF,opc12_uniqueM,
                                      opc90_rheg,opc90_uniqueF,opc90_uniqueM,
                                      bulk$deseq2$de12$genesDE,
                                      bulk$deseq2$de90$genesDE,
                                      bulk$deseq2$de150$genesDE))

umap_all_scale <- umap_all_scale[genes_DE_and_het,]
umap_all_scale <- umap_all_scale[complete.cases(umap_all_scale),]

setwd("~")
set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 2)
plot.new()
plot(embedding_umap)

umap_all.info <- c(paste0("bulk",bulk$info$day,"_",bulk$info$genotype),rep("opc12",count(opc12$info$type == "ten-cell")),rep("opc90",count(opc90$info$type == "ten-cell")))
count(opc12$info$type == "ten-cell")

pch_umap <- c(16,1,17,2,15,0,16,16)[as.numeric(factor(umap_all.info))]
col_umap <- c(rep("#000000",6),"#e41a1c","#377eb8")[as.numeric(factor(umap_all.info))]

set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 2)
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     xaxs = "i",yaxs = "i",
     xlim = c(-15,15),ylim = c(-20,15))

set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 7)
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     xaxs = "i",yaxs = "i",las = 1,
     xlim = c(-2,2),ylim = c(-6,8))

set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 12)
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     xaxs = "i",yaxs = "i",las = 1,
     xlim = c(-4,4),ylim = c(-2,2))

set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 14)
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     xaxs = "i",yaxs = "i",las = 1,
     xlim = c(-4,4),ylim = c(-2,2))











umap_all_scale <- cbind(umap_bulk_scale,umap_10c_scale)

genes_DE_and_het <- Reduce(union,list(opc12_rheg,
                                      opc90_rheg,
                                      bulk$deseq2$de12$genesDE,
                                      bulk$deseq2$de90$genesDE,
                                      bulk$deseq2$de150$genesDE))

umap_all_scale <- umap_all_scale[genes_DE_and_het,]
umap_all_scale <- umap_all_scale[complete.cases(umap_all_scale),]

setwd("~")
set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 2)


pdf(file = "./plots/umap_all_nN_2.pdf",width = 6,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,12))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,)
     # xlim = c(-4,4),ylim = c(-2,2))
legend(x = 14,y = 8,legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()





set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 7)
pdf(file = "./plots/umap_all_nN_7.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
    xlim = c(-15,5),ylim = c(-10,6))
legend(x = "topleft",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()


set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 12)
pdf(file = "./plots/umap_all_nN_12.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,5),ylim = c(-3,2))
legend(x = "topright",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()






set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 20)
pdf(file = "./plots/umap_all_nN_20.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
legend(x = "topleft",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()





set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 40)
pdf(file = "./plots/umap_all_nN_40.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
legend(x = "topleft",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()






set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 70)
pdf(file = "./plots/umap_all_nN_70.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-2,3))
legend(x = "topright",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()







set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 100)
pdf(file = "./plots/umap_all_nN_100.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-2,3))
legend(x = "topright",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()



set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 2)
pdf(file = "./plots/umap_all_nN_2.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-15,15),ylim = c(-15,15))
legend(x = "topright",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()






















set.seed(0)
embedding_umap <- uwot::umap(X = t(umap_all_scale),n_neighbors = 20)
pdf(file = "./plots/umap_all_nN_20.pdf",width = 4,height = 4,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap,
     pch = pch_umap,
     col = col_umap,
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
legend(x = "topleft",legend = c("bulk12-WT","bulk90-WT","bulk150-WT","bulk12-CKO","bulk90-CKO","bulk150-CKO","10c-12","10c-90"),xpd = T,
       pch = c(1,0,2,16,15,17,16,16),
       col = c(rep("#000000",6),"#e41a1c","#377eb8"))
dev.off()



pdf(file = "./plots/umap_all_nN_20_bulk150.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap[umap_all.info == "bulk150_WT" | umap_all.info == "bulk150_CKO",],
     pch = pch_umap[umap_all.info == "bulk150_WT" | umap_all.info == "bulk150_CKO"],
     col = col_umap[umap_all.info == "bulk150_WT" | umap_all.info == "bulk150_CKO"],
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
arrows(x0 = median(embedding_umap[umap_all.info == "bulk150_WT",1]),
       y0 = median(embedding_umap[umap_all.info == "bulk150_WT",2]),
       x1 = median(embedding_umap[umap_all.info == "bulk150_CKO",1]),
       y1 = median(embedding_umap[umap_all.info == "bulk150_CKO",2]),lwd = 2)
legend(x = "bottomright",legend = c("WT","N1P"),pch = c(2,17))
dev.off()

pdf(file = "./plots/umap_all_nN_20_bulk12.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap[umap_all.info == "bulk12_WT" | umap_all.info == "bulk12_CKO",],
     pch = pch_umap[umap_all.info == "bulk12_WT" | umap_all.info == "bulk12_CKO"],
     col = col_umap[umap_all.info == "bulk12_WT" | umap_all.info == "bulk12_CKO"],
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
arrows(x0 = median(embedding_umap[umap_all.info == "bulk12_WT",1]),
       y0 = median(embedding_umap[umap_all.info == "bulk12_WT",2]),
       x1 = median(embedding_umap[umap_all.info == "bulk12_CKO",1]),
       y1 = median(embedding_umap[umap_all.info == "bulk12_CKO",2]),lwd = 2)
legend(x = "bottomright",legend = c("WT","N1P"),pch = c(1,16))
dev.off()


pdf(file = "./plots/umap_all_nN_20_bulk90.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap[umap_all.info == "bulk90_WT" | umap_all.info == "bulk90_CKO",],
     pch = pch_umap[umap_all.info == "bulk90_WT" | umap_all.info == "bulk90_CKO"],
     col = col_umap[umap_all.info == "bulk90_WT" | umap_all.info == "bulk90_CKO"],
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
arrows(x0 = median(embedding_umap[umap_all.info == "bulk90_WT",1]),
       y0 = median(embedding_umap[umap_all.info == "bulk90_WT",2]),
       x1 = median(embedding_umap[umap_all.info == "bulk90_CKO",1]),
       y1 = median(embedding_umap[umap_all.info == "bulk90_CKO",2]),lwd = 2)
legend(x = "bottomright",legend = c("WT","N1P"),pch = c(1,16))
dev.off()


pdf(file = "./plots/umap_all_nN_20_opc.pdf",width = 2,height = 2,pointsize = 6,useDingbats = F)
par(mar = c(3.1,3.1,1,1))
plot(x = embedding_umap[umap_all.info == "opc12" | umap_all.info == "opc90",],
     pch = c(rep(1,56),rep(16,56)),
     #col = col_umap[umap_all.info == "opc12" | umap_all.info == "opc90"],
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
legend(x = "bottomright",legend = c("12 dpi","90 dpi"),pch = c(1,16))
dev.off()


# clustering bulk 90

cb90 <- embedding_umap[umap_all.info == "bulk90_CKO",]
km90 <- kmeans(x = cb90,centers = 2)
plot(x = embedding_umap[umap_all.info == "bulk90_WT",],
     pch = pch_umap[umap_all.info == "bulk90_WT"],
     col = col_umap[umap_all.info == "bulk90_WT"],
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))
points(x = cb90,
       pch = km90$cluster+15)







plot(x = embedding_umap[umap_all.info == "opc12" | umap_all.info == "opc90",],
     pch = c(rep(1,56),rep(16,56)),
     #col = col_umap[umap_all.info == "opc12" | umap_all.info == "opc90"],
     xlab = "UMAP-1",ylab = "UMAP-2",
     frame = F,
     las = 1,
     xlim = c(-3,3),ylim = c(-3,3))

pch_umap <- c(16,1,16,1,16,1,16,16)[as.numeric(factor(umap_all.info))]
{
  pdf(file = "./plots/umap_all_nN_20_bulk150_final.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F)
  par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
  plot(x = embedding_umap[umap_all.info == "bulk150_WT" | umap_all.info == "bulk150_CKO",],
       pch = pch_umap[umap_all.info == "bulk150_WT" | umap_all.info == "bulk150_CKO"],
       col = col_umap[umap_all.info == "bulk150_WT" | umap_all.info == "bulk150_CKO"],
       xlab = "UMAP-1",ylab = "UMAP-2",
       frame = F,
       las = 1,xaxs = "i",yaxs = "i",
       xlim = c(-3,3),ylim = c(-3,3),lwd = 0.5)
  lines(x = c(median(embedding_umap[umap_all.info == "bulk150_WT",1]),median(embedding_umap[umap_all.info == "bulk150_CKO",1])),
         y = c(median(embedding_umap[umap_all.info == "bulk150_WT",2]),median(embedding_umap[umap_all.info == "bulk150_CKO",2])))
  legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5)
  dev.off()
  
  pdf(file = "./plots/umap_all_nN_20_bulk12_final.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F)
  par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
  plot(x = embedding_umap[umap_all.info == "bulk12_WT" | umap_all.info == "bulk12_CKO",],
       pch = pch_umap[umap_all.info == "bulk12_WT" | umap_all.info == "bulk12_CKO"],
       col = col_umap[umap_all.info == "bulk12_WT" | umap_all.info == "bulk12_CKO"],
       xlab = "UMAP-1",ylab = "UMAP-2",
       frame = F,
       las = 1,xaxs = "i",yaxs = "i",
       xlim = c(-3,3),ylim = c(-3,3),lwd = 0.5)
  lines(x = c(median(embedding_umap[umap_all.info == "bulk12_WT",1]),median(embedding_umap[umap_all.info == "bulk12_CKO",1])),
        y = c(median(embedding_umap[umap_all.info == "bulk12_WT",2]),median(embedding_umap[umap_all.info == "bulk12_CKO",2])))
  legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5)
  dev.off()
  
  pdf(file = "./plots/umap_all_nN_20_bulk90_final.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F)
  par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
  plot(x = embedding_umap[umap_all.info == "bulk90_WT" | umap_all.info == "bulk90_CKO",],
       pch = pch_umap[umap_all.info == "bulk90_WT" | umap_all.info == "bulk90_CKO"],
       col = col_umap[umap_all.info == "bulk90_WT" | umap_all.info == "bulk90_CKO"],
       xlab = "UMAP-1",ylab = "UMAP-2",
       frame = F,
       las = 1,xaxs = "i",yaxs = "i",
       xlim = c(-3,3),ylim = c(-3,3),lwd = 0.5)
  lines(x = c(median(embedding_umap[umap_all.info == "bulk90_WT",1]),median(embedding_umap[umap_all.info == "bulk90_CKO",1][km90$cluster == 1])),
        y = c(median(embedding_umap[umap_all.info == "bulk90_WT",2]),median(embedding_umap[umap_all.info == "bulk90_CKO",2][km90$cluster == 1])))
  lines(x = c(median(embedding_umap[umap_all.info == "bulk90_WT",1]),median(embedding_umap[umap_all.info == "bulk90_CKO",1][km90$cluster == 2])),
        y = c(median(embedding_umap[umap_all.info == "bulk90_WT",2]),median(embedding_umap[umap_all.info == "bulk90_CKO",2][km90$cluster == 2])))
  legend(x = "topright",legend = c("WT","N1P"),pch = c(1,16),pt.lwd = 0.5)
  dev.off()
  
  pdf(file = "./plots/umap_all_nN_20_tencell_final.pdf",width = 1.5625,height = 1.5625,pointsize = 7,useDingbats = F)
  par(mai = c(0.25,0.25,0.05,0.05),mgp = c(1.25,0.5,0))
  plot(x = embedding_umap[umap_all.info == "opc12" | umap_all.info == "opc90",],
       pch = c(rep(1,56),rep(16,56)),
       xlab = "UMAP-1",ylab = "UMAP-2",
       frame = F,
       las = 1,xaxs = "i",yaxs = "i",
       xlim = c(-3,3),ylim = c(-3,3),lwd = 0.5)
  legend(x = "topright",legend = c("12 dpi","90 dpi"),pch = c(1,16),pt.lwd = 0.5)
  dev.off()
}
