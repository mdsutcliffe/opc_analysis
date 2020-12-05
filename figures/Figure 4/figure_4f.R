# Figure 4F

library(pheatmap)
library(RColorBrewer)
library(DESeq2)

f <- "./figures/Figure 4/CIBERSORTx_opc90.txt"
res <- read.table(file = f,header = T,sep = "\t",row.names = 1)

res <- res[,1:(which(names(res) == "P.value") - 1)]

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
              silent = T)

o_group <- unlist(as.dendrogram(p$tree_row)[[1]])
n_group <- unlist(as.dendrogram(p$tree_row)[[2]][[1]])
e_group <- c(unlist(as.dendrogram(p$tree_row)[[2]][[2]][[2]][[1]]),unlist(as.dendrogram(p$tree_row)[[2]][[2]][[2]][[2]][[1]]))
u_group <- unlist(as.dendrogram(p$tree_row)[[2]][[2]][[2]][[2]][[2]][[2]])

group <- rep("",56)
group[c(o_group,e_group,n_group)]  <- "OEN"
group[u_group] <- "U"

annotation <- data.frame(row.names = row.names(res),
                         group = group)

source("./import/import_opc90.R")

deseq_opc90 <- cbind(data.frame(row.names = opc90$rsem$symbol),opc90$rsem[,c(9+which(opc90$info$type == "ten-cell"))])
deseq_opc90 <- deseq_opc90[rowSums(deseq_opc90) > 0,]
deseq_opc90 <- deseq_opc90[,annotation$group %in% c("U","OEN")]
deseq_info <- annotation[annotation$group %in% c("U","OEN"),,drop = F]
deseq_info$group <- relevel(factor(deseq_info$group),ref = "OEN")

deseq_opc90 <- round(deseq_opc90)
deseq_opc90 <- as.matrix(deseq_opc90)
mode(deseq_opc90) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = deseq_opc90,colData = deseq_info,design = ~group)
dds <- DESeq(object = dds)

res <- results(object = dds)
resOrdered <- res[order(res$padj),]

goi <- c("Rad51c","Slx1b","Ercc4","Upf3b")

violin_groups <- c(paste(goi,"O-E-N",sep = "_"),paste(goi,"U",sep = "_"))
violin_groups <- violin_groups[c(1,5,2,6,3,7,4,8)]

expr <-apply(X = opc90$log2[match(goi,opc90$log2$symbol),((10+40):ncol(opc90$log2))[annotation$group %in% c("U","OEN")]],
             MARGIN = 1,
             FUN = function(x) c(x[deseq_info$group == "OEN"],x[deseq_info$group == "U"]))
df <- data.frame(x = rep(violin_groups,times = rep(c(sum(annotation$group == "OEN"),sum(annotation$group == "U")),4)),
                 y = c(as.numeric(c(expr[,1],expr[,2],expr[,3],expr[,4]))),
                 stringsAsFactors = F)
df$x <- factor(df$x,levels = rev(unique(df$x)))

pdf(file = "./Figure 4/figure_4f.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = df$y,y = as.numeric(df$x)+runif(n = nrow(df),min = -0.15,max = 0.15),
     pch = 16,
     cex = 0.5,
     col = "#bdbdbd",
     xlim = c(0,12),ylim = c(0.5,8.5),
     axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "Log2(TPM+1)",ylab = "",main = "Differential expression analysis",font.main = 1,cex.main = 8/7)
axis(side = 1,lwd = 0.5/0.75)
axis(side = 2,at = 0.5:8.5,labels = NA,lwd = 0.5/0.75)
axis(side = 2,at = 1:8,labels = rep(x = rev(c("O-E-N","U")),4),tick = F,las = 1,mgp = c(0,0.5,0))
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
text(x = rep(9.5,4),y = seq(1.5,7.5,2),labels = rev(goi),adj = c(0,0.5),col = rev(brewer.pal(n = 4,name = "Dark2")))
dev.off()