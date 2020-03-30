
library(pheatmap)
library(RColorBrewer)
library(DESeq2)

f <- "./Figure 4/CIBERSORTx_opc90.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)
resratio <- res[1:(which(names(res) == "P.value") - 1)] / res$Absolute.score..sig.score.

res_opc90 <- res[,1:(which(names(res) == "P.value") - 1)]

# pdf(file = "./plots/cibersort_opc90_heatmap_4E.pdf",width = 3,height = 2.25,pointsize = 7,family = "ArialMT")
p <- pheatmap(mat = res_opc90[,c(6,1,2,4,5,3)],
              # annotation_row = annotation,
              color = brewer.pal(9,"Reds"),
              clustering_method = "ward.D2",
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 50,
              treeheight_col = 10,
              show_rownames = F)

o_group <- unlist(as.dendrogram(p$tree_row)[[1]])
n_group <- unlist(as.dendrogram(p$tree_row)[[2]][[1]])
m_group <- unlist(as.dendrogram(p$tree_row)[[2]][[2]][[1]])
e_group <- c(unlist(as.dendrogram(p$tree_row)[[2]][[2]][[2]][[1]]),unlist(as.dendrogram(p$tree_row)[[2]][[2]][[2]][[2]][[1]]))
u_group <- unlist(as.dendrogram(p$tree_row)[[2]][[2]][[2]][[2]][[2]][[2]])

annotation <- data.frame(row.names = row.names(res_opc90),
                         O = (1:56) %in% o_group,
                         E = (1:56) %in% e_group,
                         N = (1:56) %in% n_group,
                         M = (1:56) %in% m_group,
                         U = (1:56) %in% u_group)
annotation_char <- data.frame(apply(X = annotation,MARGIN = 2,FUN = function(x) {
  y <- as.character(x)
  names(y) <- row.names(annotation)
  y
}))
annotation_colors <- list(O = c("TRUE" = ))

a <- data.frame(row.names = row.names(res_opc90),
                group = rep("",56),stringsAsFactors = F)
a$group[o_group] <- "O"
a$group[e_group] <- "E"
a$group[n_group] <- "N"
a$group[m_group] <- "M"
a$group[u_group] <- "U"

p <- pheatmap(mat = res_opc90,
              annotation_row = a,
              color = brewer.pal(9,"Reds"),
              clustering_method = "ward.D2",
              border_color = NA,
              lwd = 0.5/0.75,
              treeheight_row = 50,
              treeheight_col = 10,
              show_rownames = F)
a_other <- a
a_other$group[a_other$group %in% c("O","E","N")] <- "other"

source("./opc90/import_opc90.R")

deseq_opc90 <- cbind(data.frame(row.names = opc90$rsem$symbol),opc90$rsem[,c(9+which(opc90$info$type == "ten-cell"))])
deseq_opc90 <- deseq_opc90[rowSums(deseq_opc90) > 0,]
deseq_opc90 <- deseq_opc90[,a_other$group %in% c("U","other")]
deseq_info <- a_other[a_other$group %in% c("U","other"),,drop = F]
deseq_info$group <- relevel(factor(deseq_info$group),ref = "other")

deseq_opc90 <- round(deseq_opc90)
deseq_opc90 <- as.matrix(deseq_opc90)
mode(deseq_opc90) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = deseq_opc90,colData = deseq_info,design = ~group)
dds <- DESeq(object = dds)

res <- results(object = dds)
resOrdered <- res[order(res$padj),]

degenes <- row.names(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05,])
degenes_up <- row.names(resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 0,])

# goi <- c("Rad51c","Bicra","Slx1b","Upf3b","Ercc4","Nr6a1")
goi <- c("Rad51c","Slx1b","Ercc4","Upf3b")
violin_groups <- c(paste(goi,"O-E-N",sep = "_"),paste(goi,"U",sep = "_"))
violin_groups <- violin_groups[c(1,5,2,6,3,7,4,8)]

expr <-apply(X = opc90$log2[match(goi,opc90$log2$symbol),((10+40):ncol(opc90$log2))[a_other$group %in% c("U","other")]],MARGIN = 1,FUN = function(x) c(x[deseq_info$group == "other"],x[deseq_info$group == "U"]))
df <- data.frame(x = rep(violin_groups,times = rep(c(sum(a_other$group == "other"),sum(a_other$group == "U")),4)),
                 y = c(as.numeric(c(expr[,1],expr[,2],expr[,3],expr[,4]))),
                 stringsAsFactors = F)
df$x <- factor(df$x,levels = rev(unique(df$x)))

pdf(file = "./plots/figure_4f.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0))
plot(x = df$y,y = as.numeric(df$x)+runif(n = nrow(df),min = -0.15,max = 0.15),
     pch = 16,
     cex = 0.5,
     col = "#bdbdbd",
     xlim = c(0,12),
     ylim = c(0.5,8.5),
     axes = F,
     xaxs = "i",
     yaxs = "i",
     xlab = "Log2(TPM+1)",ylab = "")
vioplot::vioplot(formula = y ~ x,
                 data = df,
                 drawRect = F,
                 horizontal = T,
                 col = "#00000000",
                 border = rev(rep(x = brewer.pal(n = 4,name = "Dark2"),each = 2)),
                 names = rep(x = rev(c("O-E-N","U")),4),
                 lwd = 0.5/0.75,
                 axes = F,
                 add = T)
sapply(1:length(levels(df$x)),function(i) {
  m <- mean(df$y[df$x == levels(df$x)[i]])
  lines(x = rep(m,2),y = c(i-0.25,i+0.25),lwd = 1/0.75)
  s <- sd(df$y[df$x == levels(df$x)[i]])/sqrt(length(df$y[df$x == levels(df$x)[i]]))
  lines(x = c(m-s,m+s),y = rep(i,2),lwd = 0.5/0.75)
  lines(x = rep(m-s,2),y = c(i-0.15,i+0.15),lwd = 0.5/0.75)
  lines(x = rep(m+s,2),y = c(i-0.15,i+0.15),lwd = 0.5/0.75)
})
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = 0.5:8.5,labels = NA,lwd = 0.5/0.75)
axis(side = 2,at = 1:8,labels = rep(x = rev(c("O-E-N","U")),4),tick = F,las = 1,mgp = c(0,0.5,0))
text(x = rep(9.5,4),y = seq(1.5,7.5,2),labels = rev(goi),adj = c(0,0.5),col = rev(brewer.pal(n = 4,name = "Dark2")))
dev.off()

pheatmap(opc90)
expr <-apply(X = opc90$log2[match(degenes,opc90$log2$symbol),((10+40):ncol(opc90$log2))[a_other$group %in% c("U","other")]],MARGIN = 1,FUN = function(x) c(x[deseq_info$group == "other"],x[deseq_info$group == "U"]))
df <- data.frame(row.names = names(opc90$info$log2[,9+which(opc90$info$type == "ten-cell")])= rep(violin_groups,times = rep(c(sum(a_other$group == "other"),sum(a_other$group == "U")),4)),
                 y = c(as.numeric(c(expr[,1],expr[,2],expr[,3],expr[,4]))),
                 stringsAsFactors = F)
df$x <- factor(df$x,levels = rev(unique(df$x)))

pdf(file = "~/ugroup_de.pdf",width = 6,height = 7)
pheatmap(mat = cbind(data.frame(row.names = opc90$log2$symbol[opc90$log2$symbol %in% degenes]),opc90$log2[opc90$log2$symbol %in% degenes,9+40+which(a_other$group %in% c("U","other"))]),
         color = rev(brewer.pal(n = 11,name = "RdBu")),breaks = seq(-4,4,length.out = 12),
         scale = "row",show_rownames = F,show_colnames = F,
         clustering_method = "ward.D2",
         annotation = data.frame(row.names = row.names(a_other)[a_other$group %in% c("U","other")],
                                 group = a_other$group[a_other$group %in% c("U","other")]))
dev.off()