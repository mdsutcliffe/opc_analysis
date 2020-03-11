# opc12 ranked dispersion

source("./opc12/RHEGs_opc12.R")

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

pdf(file = "./plots/ranked_dispersion_opc12.pdf",width = 2.25,height = 2.5,pointsize = 7,useDingbats = F)
par(mai = c(0.15,0.5,0,0),mfrow = c(2,1),mgp = c(0.85,0.6,0))
plot(x = c(1:10585,c(1:10585)[!(row.names(avg_vars_M_sort) %in% opc12_rheg)]),
     y = c(avg_vars_M_sort$pooled_var,avg_vars_M_sort$tencell_var[!(row.names(avg_vars_M_sort) %in% opc12_rheg)]),
     col = c(rep("#00000020",10585),rep("#bcbddc",count(!(row.names(avg_vars_M_sort) %in% opc12_rheg)))),
     frame = F,las = 1,axes = F,
     xlim = c(1,10585),ylim = c(0,12),
     xlab = NA,
     ylab = NA,
     lwd = 0.5)
points(x = c(1:10585)[row.names(avg_vars_M_sort) %in% opc12_rheg],y = avg_vars_M_sort$tencell_var[row.names(avg_vars_M_sort) %in% opc12_rheg],
       col = "#756bb1",
       lwd = 0.5)
axis(side = 2,las = 1,lwd = 0.5)
text(x = 10585,y = 12,labels = "Male",adj = c(1,1))
plot(x = c(1:10585,c(1:10585)[!(row.names(avg_vars_F_sort) %in% opc12_rheg)]),
     y = c(avg_vars_F_sort$pooled_var,avg_vars_F_sort$tencell_var[!(row.names(avg_vars_F_sort) %in% opc12_rheg)]),
     col = c(rep("#00000020",10585),rep("#99d8c9",count(!(row.names(avg_vars_F_sort) %in% opc12_rheg)))),
     frame = F,las = 1,axes = F,
     xlim = c(1,10585),ylim = c(0,12),
     xlab = "Male/Female gene rank",
     ylab = NA,
     lwd = 0.5)
points(x = c(1:10585)[row.names(avg_vars_F_sort) %in% opc12_rheg],y = avg_vars_F_sort$tencell_var[row.names(avg_vars_F_sort) %in% opc12_rheg],
       col = "#31a354",
       lwd = 0.5)
axis(side = 1,at = c(1,10585),lwd = 0.5)
axis(side = 2,las = 1,lwd = 0.5)
text(x = 10585,y = 12,labels = "Female",adj = c(1,1))
title(ylab = "Ranked dispersion",mgp = c(1.6,0,0))
dev.off()