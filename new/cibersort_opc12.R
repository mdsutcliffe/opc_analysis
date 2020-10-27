# Figure 4E

library(pheatmap)
library(RColorBrewer)

f <- "./new/CIBERSORTx_opc12.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)

res <- res[,1:(which(names(res) == "P.value") - 1)]
names(res)[names(res) == "MO"] <- "Myelinating\noligodendrocyte"

p <- pheatmap(mat = res,
              color = brewer.pal(9,"Reds"),
              clustering_method = "ward.D2",
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 20,
              treeheight_col = 10,
              fontsize = 7,
              show_rownames = F,
              angle_col = 90,
              filename =  "./plots/heatmap_cibersort_opc12.pdf",
              width = 3.125,
              height = 3.125)