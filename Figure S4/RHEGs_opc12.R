# Get opc12 candidate genes

source("./functions/normalizeTPM.R")
source("./opc12/import_opc12.R")

resF <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
resM <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))

nAppearancesF <- table(unlist(x = resF))
nAppearancesM <- table(unlist(x = resM))

opc12_uniqueF <- setdiff(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])
opc12_uniqueM <- setdiff(x = names(nAppearancesM)[nAppearancesM >= 75],y = names(nAppearancesF)[nAppearancesF >= 75])
opc12_rheg <- intersect(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])

opc12_rheg_mat <- opc12$log2[match(x = opc12_rheg,table = opc12$log2$symbol),c(1:9,9+which(opc12$info$type == "ten-cell"))]
row.names(opc12_rheg_mat) <- opc12_rheg_mat$symbol
opc12_rheg_mat <- opc12_rheg_mat[,10:ncol(opc12_rheg_mat)]
opc12_rheg_mat <- as.matrix(opc12_rheg_mat)

pdf(file = "./plots/figure_s4a.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.4,0,0),mgp = c(2.1,0.6,0))
hist(x = as.numeric(nAppearancesF),
     breaks = 0:100,
     col = "#bdbdbd",
     border = F,
     xaxs = "i",
     yaxs = "i",
     ylim = c(0,400),
     las = 1,
     main = NA,xlab = NA,
     lwd = 0.5/0.75)
title(xlab = "Percent of simulations",mgp = c(1.6,0.6,0))
lines(x = c(75,75,85),y = c(0,100,100),lwd = 0.5/0.75,col = "#de2d26")
text(x = 75,y = 140,labels = "Robust\nheterogeneities",col = "#de2d26")
text(x = 76,y = 75,labels = length(c(opc12_rheg,opc12_uniqueF)),adj = c(0,0.5),col = "#de2d26")
dev.off()

pdf(file = "./plots/figure_s4b.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.4,0,0),mgp = c(2.1,0.6,0))
hist(x = as.numeric(nAppearancesM),
     breaks = 0:100,
     col = "#bdbdbd",
     border = F,
     xaxs = "i",
     yaxs = "i",
     ylim = c(0,400),
     las = 1,
     lwd = 0.5/0.75,
     main = NA,xlab = NA)
title(xlab = "Percent of simulations",mgp = c(1.6,0.6,0))
lines(x = c(75,75,85),y = c(0,100,100),lwd = 0.5/0.75,col = "#de2d26")
text(x = 75,y = 140,labels = "Robust\nheterogeneities",col = "#de2d26")
text(x = 76,y = 75,labels = length(c(opc12_rheg,opc12_uniqueM)),adj = c(0,0.5),col = "#de2d26")
dev.off()