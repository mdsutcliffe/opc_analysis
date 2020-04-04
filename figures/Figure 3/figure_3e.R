# Figure 3E

fList <- list.files(path = "./data/opc12_rdatafiles",full.names = T)

arvListF <- list()
arvListM <- list()
for (i in 1:length(fList)) {
  load(fList[i])
  arvListF[[i]] <- resFemale$arv
  arvListM[[i]] <- resMale$arv
}

gene_set_F <- table(unlist(sapply(1:100,function(x) arvListF[[x]]$symbol)))
gene_set_F <- names(gene_set_F)[gene_set_F == 100]

gene_set_M <- table(unlist(sapply(1:100,function(x) arvListM[[x]]$symbol)))
gene_set_M <- names(gene_set_M)[gene_set_M == 100]

gene_set <- intersect(gene_set_F,gene_set_M)

tencell_var_F <- cbind(data.frame(row.names = gene_set),do.call(cbind,lapply(X = arvListF,FUN = function(x) x$tencell_var[x$symbol %in% gene_set])))
tencell_var_M <- cbind(data.frame(row.names = gene_set),do.call(cbind,lapply(X = arvListM,FUN = function(x) x$tencell_var[x$symbol %in% gene_set])))

pooled_var_F <- cbind(data.frame(row.names = gene_set),do.call(cbind,lapply(X = arvListF,FUN = function(x) x$pooled_var[x$symbol %in% gene_set])))
pooled_var_M <- cbind(data.frame(row.names = gene_set),do.call(cbind,lapply(X = arvListM,FUN = function(x) x$pooled_var[x$symbol %in% gene_set])))

avg_vars_F <- data.frame(cbind(tencell_var = apply(X = tencell_var_F,MARGIN = 1,FUN = mean),pooled_var = apply(X = pooled_var_F,MARGIN = 1,FUN = mean)))
avg_vars_M <- data.frame(cbind(tencell_var = apply(X = tencell_var_M,MARGIN = 1,FUN = mean),pooled_var = apply(X = pooled_var_M,MARGIN = 1,FUN = mean)))

avg_vars_F_sort <- avg_vars_F[order(avg_vars_F$tencell_var,decreasing = T),]
avg_vars_M_sort <- avg_vars_M[order(avg_vars_M$tencell_var,decreasing = T),]

avg_vars_F_sort <- log10(avg_vars_F_sort)
avg_vars_M_sort <- log10(avg_vars_M_sort)

source("./functions/opc12_rheg.R")

goi <- c("Aldh1a1","Axl","Ltbp4","Bmp1","Nbl1","Wnt7a","Ctnnb1")

pdf(file = "./Figure 3/figure_3e.pdf",width = 3.125,height = 3.125,pointsize = 7,useDingbats = F,bg = "white")
par(mfrow = c(2,1),xpd = T)

par(mai = c(0.0625,0.5,0.25,0.5),mgp = c(1.6,0.6,0))
plot(x = NA,y = NA,
     xlim = c(1,10585),ylim = c(-0.5,1),
     frame = F,las = 1,axes = F,
     xlab = NA,ylab = NA,
     xaxs = "i",yaxs = "i",
     lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = 0.5/0.75,at = seq(-0.5,1,0.5))
points(x = 1:10585,
       y = avg_vars_M_sort$pooled_var,
       col = "#00000010",pch = 16,cex = 0.5)
points(x = (1:10585)[!(row.names(avg_vars_M_sort) %in% c(opc12$candidates$RHEG,opc12$candidates$Mspecific))],
       y = avg_vars_M_sort$tencell_var[!(row.names(avg_vars_M_sort) %in% c(opc12$candidates$RHEG,opc12$candidates$Mspecific))],
       col = "#d01c8b",pch = 16,cex = 0.5)
points(x = c(1:10585)[row.names(avg_vars_M_sort) %in% setdiff(c(opc12$candidates$RHEG,opc12$candidates$Mspecific),goi)],
       y = avg_vars_M_sort$tencell_var[row.names(avg_vars_M_sort) %in% setdiff(c(opc12$candidates$RHEG,opc12$candidates$Mspecific),goi)],
       col = "#ffbf00",pch = 16,cex = 0.5)
points(x = c(1:10585)[row.names(avg_vars_M_sort) %in% goi],y = avg_vars_M_sort$tencell_var[row.names(avg_vars_M_sort) %in% goi],
       col = "#000000",pch = 16,cex = 0.5)
text(x = 10585,y = 0.5,labels = "Male",adj = c(1,1))

par(mai = c(0.25,0.5,0.0625,0.5),mgp = c(1.1,0.6,0))
plot(x = NA,y = NA,
     xlim = c(1,10585),ylim = c(-0.5,1),
     frame = F,las = 1,axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "Male/Female gene rank",ylab = NA,
     lwd = 0.5/0.75)
axis(side = 1,at = c(1,10585),lwd = 0.5/0.75)
axis(side = 2,las = 1,lwd = 0.5/0.75,at = seq(-0.5,1,0.5))
points(x = 1:10585,
       y = avg_vars_F_sort$pooled_var,
       col = "#00000010",pch = 16,cex = 0.5)
points(x = (1:10585)[!(row.names(avg_vars_F_sort) %in% c(opc12$candidates$RHEG,opc12$candidates$Fspecific))],
       y = avg_vars_F_sort$tencell_var[!(row.names(avg_vars_F_sort) %in% c(opc12$candidates$RHEG,opc12$candidates$Fspecific))],
       col = "#4dac26",pch = 16,cex = 0.5)
points(x = (1:10585)[row.names(avg_vars_F_sort) %in% setdiff(c(opc12$candidates$RHEG,opc12$candidates$Fspecific),goi)],
       y = avg_vars_F_sort$tencell_var[row.names(avg_vars_F_sort) %in% setdiff(c(opc12$candidates$RHEG,opc12$candidates$Fspecific),goi)],
       col = "#ffbf00",pch = 16,cex = 0.5)
points(x = (1:10585)[row.names(avg_vars_F_sort) %in% goi],
       y = avg_vars_F_sort$tencell_var[row.names(avg_vars_F_sort) %in% goi],
       col = "#000000",pch = 16,cex = 0.5)
text(x = 10585,y = 0.5,labels = "Female",adj = c(1,1))
title(ylab = "Ranked dispersion",mgp = c(2.4,0,0))
dev.off()