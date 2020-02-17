# Heterogeneity analysis in opc12

# File paths
opc12.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data"
opc12F.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/resample_results_opc12_female"
opc12M.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/resample_results_opc12_male"

getCandidates <- function(path,i) {
  f <- read.table(file = paste0(path,"/resample_",sprintf("%03d",i),"_variance.tsv"),sep = "\t",header = T,row.names = 1)
  genes <- row.names(f)[f$candidate]
}

opc12F_candidates <- lapply(1:100,function(i) getCandidates(path = opc12F.path,i = i))
opc12M_candidates <- lapply(1:100,function(i) getCandidates(path = opc12M.path,i = i))

opc12F_candidates_table <- table(unlist(opc12F_candidates))
opc12M_candidates_table <- table(unlist(opc12M_candidates))

opc12F_candidates_table_filter <- names(opc12F_candidates_table)[opc12F_candidates_table >= 75]
opc12M_candidates_table_filter <- names(opc12M_candidates_table)[opc12M_candidates_table >= 75]
opc12_sexIndependentCandidates <- intersect(opc12F_candidates_table_filter,opc12M_candidates_table_filter)



opc12 <- read.csv(file = paste0(opc12.path,"/rsem_opc12.csv"),stringsAsFactors = F)

opc12F_samples_var <- data.frame(symbol = opc12$symbol,candidate = opc12$symbol %in% opc12F_candidates_table_filter)
opc12M_samples_var <- data.frame(symbol = opc12$symbol,candidate = opc12$symbol %in% opc12M_candidates_table_filter)

opc12F_controls_var <- data.frame(symbol = opc12$symbol,candidate = opc12$symbol %in% opc12F_candidates_table_filter)
opc12M_controls_var <- data.frame(symbol = opc12$symbol,candidate = opc12$symbol %in% opc12M_candidates_table_filter)

getCandidateVars <- function(path,i) {
  f <- read.table(file = paste0(path,"/resample_",sprintf("%03d",i),"_variance.tsv"),sep = "\t",header = T)
  names(f) <- paste0(names(f),".",sprintf("%03d",i))
  names(f)[1] <- "symbol"
  return(f)
}

for (i in 1:100) {
  opc12F_samples_var <- merge(opc12F_samples_var,getCandidateVars(path = opc12F.path,i = i)[,c(1,2),drop = F],by = "symbol",all = T)
  opc12M_samples_var <- merge(opc12M_samples_var,getCandidateVars(path = opc12M.path,i = i)[,c(1,2),drop = F],by = "symbol",all = T)
  
  opc12F_controls_var <- merge(opc12F_controls_var,getCandidateVars(path = opc12F.path,i = i)[,c(1,3),drop = F],by = "symbol",all = T)
  opc12M_controls_var <- merge(opc12M_controls_var,getCandidateVars(path = opc12M.path,i = i)[,c(1,3),drop = F],by = "symbol",all = T)
}

sum(complete.cases(opc12F_samples_var[,3:102]))
sum(complete.cases(opc12M_samples_var[,3:102]))

opc12F_samples_genes <- opc12F_samples_var$symbol[complete.cases(opc12F_samples_var[,3:102])]
opc12M_samples_genes <- opc12M_samples_var$symbol[complete.cases(opc12M_samples_var[,3:102])]

opc12F_controls_genes <- opc12F_controls_var$symbol[complete.cases(opc12F_controls_var[,3:102])]
opc12M_controls_genes <- opc12M_controls_var$symbol[complete.cases(opc12M_controls_var[,3:102])]

geneList <- intersect(intersect(opc12F_samples_genes,opc12M_samples_genes),
                      intersect(opc12F_controls_genes,opc12M_controls_genes))

opc12F_samples_var$var_mean <- apply(opc12F_samples_var[,3:102],1,mean)
opc12M_samples_var$var_mean <- apply(opc12M_samples_var[,3:102],1,mean)

opc12F_controls_var$var_mean <- apply(opc12F_controls_var[,3:102],1,mean)
opc12M_controls_var$var_mean <- apply(opc12M_controls_var[,3:102],1,mean)

opc12F_samples_var_genes <- opc12F_samples_var[opc12F_samples_var$symbol %in% geneList,]
opc12M_samples_var_genes <- opc12M_samples_var[opc12M_samples_var$symbol %in% geneList,]

opc12F_controls_var_genes <- opc12F_controls_var[opc12F_controls_var$symbol %in% geneList,]
opc12M_controls_var_genes <- opc12M_controls_var[opc12M_controls_var$symbol %in% geneList,]

rank_opc12F_samples <- order(opc12F_samples_var_genes$var_mean)
rank_opc12M_samples <- order(opc12M_samples_var_genes$var_mean)

pdf(file = "plots/opc12F_rankedDispersion.pdf",width = 8,height =8)
plot(1:9139,rev(opc12F_controls_var_genes$var_mean[rank_opc12F_samples]),
     pch = 16,col = "#44444444",
     xlab = "Gene rank",ylab = "Mean dispersion",axes = F, xlim = c(0,9140),ylim = c(0,6),
     main = "Female ranked dispersion")
points(1:9139,rev(opc12F_samples_var_genes$var_mean[rank_opc12F_samples]),
       pch = 16,col = "#de2d26")
axis(1,at = c(1,9139))
axis(2,at = 0:6,las = 1)
legend(x = "topright",legend = c("split-pool controls","ten-cell samples"),col = c("#44444444","#de2d26"),pch=c(16,16),)
dev.off()

pdf(file = "plots/opc12M_rankedDispersion.pdf",width = 8,height =8)
plot(1:9139,rev(opc12M_controls_var_genes$var_mean[rank_opc12M_samples]),
     pch = 16,col = "#44444444",
     xlab = "Gene rank",ylab = "Mean dispersion",axes = F, xlim = c(0,9140),ylim = c(0,6),
     main = "Male ranked dispersion")
points(1:9139,rev(opc12M_samples_var_genes$var_mean[rank_opc12M_samples]),
       pch = 16,col = "#de2d26")
axis(1,at = c(1,9139))
axis(2,at = 0:6,las = 1)
legend(x = "topright",legend = c("split-pool controls","ten-cell samples"),col = c("#44444444","#de2d26"),pch=c(16,16),)
dev.off()

ratioMF_samples <- opc12F_samples_var_genes$var_mean / opc12M_samples_var_genes$var_mean
ratioMF_controls <- opc12F_controls_var_genes$var_mean / opc12M_controls_var_genes$var_mean
rank_ratioMF_samples <- order(ratioMF_samples)

# plot(1:9139,rev(ratioMF_controls[rank_ratioMF_samples]),
#      pch = 16,col = "#44444444",
#      xlab = "Gene rank",ylab = "Female/Male mean dispersion",axes = F, xlim = c(0,9140),ylim = c(0,4),
#      main = "Female/Male ranked dispersion")
# points(1:9139,rev(ratioMF_samples[rank_ratioMF_samples]),
#        pch = 16,col = "#de2d2644")
# x <- 1:9139
# y = rev(ratioMF_samples[rank_ratioMF_samples])
# points(x[geneList %in% opc12_sexIndependentCandidates],y[geneList %in% opc12_sexIndependentCandidates],pch = 16,col = "#3182bd")
# axis(1,at = c(1,9139))
# axis(2,at = 0:4,las = 1)

pdf(file = "plots/opc12F_rankedDispersion_annotated.pdf",width = 8,height =8)
plot(1:9139,rev(opc12F_controls_var_genes$var_mean[rank_opc12F_samples]),
     pch = 16,col = "#44444422",
     xlab = "Gene rank",ylab = "Mean dispersion",axes = F, xlim = c(0,9140),ylim = c(0,6),
     main = "Female ranked dispersion")
points(1:9139,rev(opc12F_samples_var_genes$var_mean[rank_opc12F_samples]),
       pch = 16,col = "#fc927222")
points(which(rev(opc12F_samples_var_genes$var_mean[rank_opc12F_samples]) %in% opc12F_samples_var_genes$var_mean[opc12F_samples_var_genes$symbol %in% opc12_sexIndependentCandidates]),
       rev(sort(opc12F_samples_var_genes$var_mean[opc12F_samples_var_genes$symbol %in% opc12_sexIndependentCandidates])),
       pch = 21,col = "#de2d26")
axis(1,at = c(1,9139))
axis(2,at = 0:6,las = 1)
legend(x = "topright",legend = c("split-pool controls","ten-cell samples","RHEGs"),col = c("#44444444","#fc927244","#de2d26"),pch=c(16,16,21),)
dev.off()

pdf(file = "plots/opc12M_rankedDispersion_annotated.pdf",width = 8,height =8)
plot(1:9139,rev(opc12M_controls_var_genes$var_mean[rank_opc12M_samples]),
     pch = 16,col = "#44444422",
     xlab = "Gene rank",ylab = "Mean dispersion",axes = F, xlim = c(0,9140),ylim = c(0,6),
     main = "Male ranked dispersion")
points(1:9139,rev(opc12M_samples_var_genes$var_mean[rank_opc12M_samples]),
       pch = 16,col = "#fc927222")
points(which(rev(opc12M_samples_var_genes$var_mean[rank_opc12M_samples]) %in% opc12M_samples_var_genes$var_mean[opc12M_samples_var_genes$symbol %in% opc12_sexIndependentCandidates]),
       rev(sort(opc12M_samples_var_genes$var_mean[opc12M_samples_var_genes$symbol %in% opc12_sexIndependentCandidates])),
       pch = 21,col = "#de2d26")
axis(1,at = c(1,9139))
axis(2,at = 0:6,las = 1)
legend(x = "topright",legend = c("split-pool controls","ten-cell samples","RHEGs"),col = c("#44444444","#fc927244","#de2d26"),pch=c(16,16,21),)
dev.off()