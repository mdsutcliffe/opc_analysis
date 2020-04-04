# Supplementary Figure S4

source("./functions/opc12_rheg.R")

resF <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
resM <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))

nAppearancesF <- table(unlist(x = resF))
nAppearancesM <- table(unlist(x = resM))

pdf(file = "./Figure S4/figure_s4a.pdf",width = 3.125,height = 3.125,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(2.1,0.6,0),lwd = 0.5/0.75)
hist(x = as.numeric(nAppearancesF),
     breaks = 0:100,
     col = "#636363",
     border = NA,
     xaxs = "i",
     yaxs = "i",
     ylim = c(0,400),
     las = 1,
     main = NA,xlab = NA,
     lwd = 0.5/0.75)
title(xlab = "Percent of simulations",mgp = c(1.6,0.6,0))
lines(x = c(75,75,85),y = c(0,100,100),lwd = 0.5/0.75,col = "#de2d26")
text(x = 75,y = 140,labels = "Robust\nheterogeneities",col = "#de2d26")
text(x = 76,y = 75,labels = length(c(opc12$candidates$RHEG,opc12$candidates$Fspecific)),adj = c(0,0.5),col = "#de2d26")
dev.off()

pdf(file = "./Figure S4/figure_s4b.pdf",width = 3.125,height = 3.125,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(2.1,0.6,0),lwd = 0.5/0.75)
hist(x = as.numeric(nAppearancesM),
     breaks = 0:100,
     col = "#636363",
     border = NA,
     xaxs = "i",
     yaxs = "i",
     ylim = c(0,400),
     las = 1,
     lwd = 0.5/0.75,
     main = NA,xlab = NA)
title(xlab = "Percent of simulations",mgp = c(1.6,0.6,0))
lines(x = c(75,75,85),y = c(0,100,100),lwd = 0.5/0.75,col = "#de2d26")
text(x = 75,y = 140,labels = "Robust\nheterogeneities",col = "#de2d26")
text(x = 76,y = 75,labels = length(c(opc12$candidates$RHEG,opc12$candidates$Mspecific)),adj = c(0,0.5),col = "#de2d26")
dev.off()