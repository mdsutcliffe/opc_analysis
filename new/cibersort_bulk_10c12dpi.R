f.cibersort_bulk <- "./figures/Figure 2/CIBERSORTx_opcBulk.txt"
f.cibersort_opc12 <- "./new/CIBERSORTx_opc12.txt"

cibersort_bulk <- read.table(file = f.cibersort_bulk,header = T,sep = "\t",row.names = 1)
cibersort_opc12 <- read.table(file = f.cibersort_opc12,header = T,sep = "\t",row.names = 1)

cibersort_bulk <- cibersort_bulk[,1:(which(names(cibersort_bulk) == "P.value") - 1)]
cibersort_opc12 <- cibersort_opc12[,1:(which(names(cibersort_opc12) == "P.value") - 1)]

names(cibersort_bulk)[names(cibersort_bulk) == "MO"] <- "Myelinating\noligodendrocyte"
names(cibersort_opc12)[names(cibersort_opc12) == "MO"] <- "Myelinating\noligodendrocyte"

cs_combined <- rbind(cibersort_bulk,cibersort_opc12)

library(pheatmap)
library(RColorBrewer)

p <- pheatmap(mat = cibersort_opc12[,c(2,6,4,5,1,3)],
              color = brewer.pal(9,"Reds"),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./new/cibersort_revisions_opc12.pdf",
              width = 3.125,
              height = 3.125)


wt_12dpi <- cibersort_bulk[grepl(pattern = "12dpi_WT",x = row.names(cibersort_bulk)),]
cko_12dpi <- cibersort_bulk[grepl(pattern = "12dpi_CKO",x = row.names(cibersort_bulk)),]
wt_90dpi <- cibersort_bulk[grepl(pattern = "90dpi_WT",x = row.names(cibersort_bulk)),]
cko_90dpi <- cibersort_bulk[grepl(pattern = "90dpi_CKO",x = row.names(cibersort_bulk)),]
wt_150dpi <- cibersort_bulk[grepl(pattern = "150dpi_WT",x = row.names(cibersort_bulk)),]
cko_150dpi <- cibersort_bulk[grepl(pattern = "150dpi_Tumor",x = row.names(cibersort_bulk)),]

p <- pheatmap(mat = wt_12dpi[,c(2,6,4,5,1,3)],
              color = brewer.pal(9,"Reds"),breaks = seq(0,12,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./new/cibersort_revisions_bulk_wt12.pdf",
              width = 2.25,
              height = 1.5625)

p <- pheatmap(mat = cko_12dpi[,c(2,6,4,5,1,3)],
              color = brewer.pal(9,"Reds"),breaks = seq(0,12,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./new/cibersort_revisions_bulk_cko12.pdf",
              width = 2.25,
              height = 1.5625)

p <- pheatmap(mat = wt_90dpi[,c(2,6,4,5,1,3)],
              color = brewer.pal(9,"Reds"),breaks = seq(0,12,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./new/cibersort_revisions_bulk_wt90.pdf",
              width = 2.25,
              height = 1.5625)

p <- pheatmap(mat = cko_90dpi[,c(2,6,4,5,1,3)],
              color = brewer.pal(9,"Reds"),breaks = seq(0,12,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./new/cibersort_revisions_bulk_cko90.pdf",
              width = 2.25,
              height = 1.5625)


p <- pheatmap(mat = wt_150dpi[,c(2,6,4,5,1,3)],
              color = brewer.pal(9,"Reds"),breaks = seq(0,12,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./new/cibersort_revisions_bulk_wt150.pdf",
              width = 2.25,
              height = 1.5625)

p <- pheatmap(mat = cko_150dpi[,c(2,6,4,5,1,3)],
              color = brewer.pal(9,"Reds"),breaks = seq(0,12,length.out = 10),
              clustering_method = "ward.D2",cluster_cols = F,
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./new/cibersort_revisions_bulk_cko150.pdf",
              width = 2.25,
              height = 1.5625)