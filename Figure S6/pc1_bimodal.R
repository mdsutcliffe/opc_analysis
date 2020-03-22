# PCA subtype analysis from PMID 28697342

source("./opc12/import_opc12.R")

load("./build/orthologs.RData")

pc <- read.csv("./external/219026_2_supp_5782669_py1bdv.csv",stringsAsFactors = F)[,1:10]

opc12_pc <- cbind(opc12$tpm[,1:9],log2(opc12$tpm[,10:ncol(opc12$tpm)] / 100 + 1))

opc12_pc$symbol <- orthologs$HGNC.symbol[match(x = opc12_pc$symbol,table = orthologs$MGI.symbol)]

opc12_pc <- opc12_pc[match(x = pc$Gene,table = opc12_pc$symbol),]

opc12_pc1 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC1 * x)
opc12_pc2 <- apply(X = opc12_pc[,10:ncol(opc12_pc)],MARGIN = 2,FUN = function(x) pc$PC2 * (x - (pc$PC1 * x)))

opc12_projection1 <- colSums(x = opc12_pc1,na.rm = T)
opc12_projection2 <- colSums(x = opc12_pc2,na.rm = T)


# THIS ONE
pdf(file = "./plots/figure_s6a.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
hist(x = opc12_projection1[opc12$info$type == "ten-cell"],
     breaks = seq(-30,30,length.out = 23),
     col = "#969696",
     border = NA,
     xlim = c(-30,30),ylim = c(0,8),
     xaxs = "i",yaxs = "i",
     xlab = "PC1",ylab = "Frequency",main = NA,
     axes = F)
axis(side = 1,at = seq(-30,30,10),labels = c(-30,NA,NA,0,NA,NA,30),lwd = 0.5/0.75,mgp = c(1.6,0.6,0))
axis(side = 2,at = 0:8,lwd = 0.5/0.75,mgp = c(1.6,0.6,0),las = 1)
lines(x = c(-2.55,-26),y = rep(-0.6,2),lwd = 0.5/0.75)
lines(x = c(2.5,26),y = rep(-0.6,2),lwd = 0.5/0.75)
text(x = -27.5/2,y = -0.8571429,labels = "Proneural",cex = 6/7)
text(x = 27.5/2,y = -0.8571429,labels = "Mesenchymal",cex = 6/7)
dev.off()

# expanded x range
pdf(file = "./plots/figure_s6a.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
hist(x = opc12_projection1[opc12$info$type == "ten-cell"],
     breaks = seq(-30,30,length.out = 23),
     col = "#969696",
     border = NA,
     xlim = c(-30,30)*2,ylim = c(0,8),
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Frequency",main = NA,
     axes = F)
axis(side = 1,at = seq(-30,30,10)*2,labels = c(-30,NA,NA,0,NA,NA,30)*2,lwd = 0.5/0.75,mgp = c(1.6,0.6,0))
axis(side = 2,at = 0:8,lwd = 0.5/0.75,mgp = c(1.6,0.6,0),las = 1)
lines(x = c(-2.55,-26)*2,y = rep(-0.6,2),lwd = 0.5/0.75)
lines(x = c(2.5,26)*2,y = rep(-0.6,2),lwd = 0.5/0.75)
text(x = -27.5,y = -0.8571429,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -0.8571429,labels = "Mesenchymal",cex = 6/7)
dev.off()

shift <- 120
x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"]) + shift

sapply(list.files("./external/Stochprof/R",full.names = T),source)
library(numDeriv)
library(MASS)

ml <- stochasticProfilingML()

# 3 x 1 no 

# Maximum likelihood estimate (MLE):
#     p_1  mu_1_gene_Gene 1  mu_2_gene_Gene 1     sigma 
#  0.8453            2.5930            1.5360    0.0790 
# 
# Value of negative log-likelihood function at MLE:
#   212.2336 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   440.5686 
# 
# Approx. 95% confidence intervals for MLE:
#                      lower     upper
# p_1              0.8018084 0.8806681
# mu_1_gene_Gene 1 2.5734546 2.6125454
# mu_2_gene_Gene 1 1.3432111 1.7287889
# sigma            0.0548461 0.1137911
sapply(X = 1:10,FUN = function(i) {
  plot(0:100,dnorm(x = 0:100,ml$mle[2],ml$mle[4]^2,log = T) + dnorm(x = 0:100,ml$mle[3],ml$mle[4]^2,log = T))
})
d.sum.of.mixtures.LNLN()

xvals <- seq(-60,60,0.1)
y1 <- dnorm(x = xvals,mean = 10^(ml$mle[3] - log10(shift)),sd = 10^(ml$mle[4])) * (1 - ml$mle[1])
y2 <- dnorm(x = xvals,mean = 10^(ml$mle[2] - log10(shift)),sd = 10^(ml$mle[4])) * ml$mle[1]
plot(x = NA,y = NA,
     xlim = range(xvals),ylim = c(0,1))
lines(x = xvals,y = y1)
lines(x = xvals,y = y2)






# Agglomerations
set.seed(0)
nCells <- sample(x = seq(700,1000,10),size = 100,replace = T)
agglomerate_pc <- t(sapply(X = nCells,FUN = function(nC) {
  cellGroups <- sample(x = 9 + which(opc12$info$type == "ten-cell"),size = nC/10,replace = T)
  
  iOPC12 <- cbind(opc12$tpm[,1:9],log2(rowMeans(opc12$tpm[,cellGroups]) / 100 + 1))
  
  iOPC12$symbol <- orthologs$HGNC.symbol[match(x = iOPC12$symbol,table = orthologs$MGI.symbol)]
  
  iOPC12 <- iOPC12[match(x = pc$Gene,table = iOPC12$symbol),]
  
  iOPC12 <- iOPC12[,10:ncol(iOPC12)]
  
  iPC1 <- pc$PC1 * iOPC12
  iPC2 <- pc$PC2 * (iOPC12 - (pc$PC1 * iOPC12))
  
  iProjection1 <- sum(x = iPC1,na.rm = T)
  iProjection2 <- sum(x = iPC2,na.rm = T)
  
  return(c(iProjection1,iProjection2))
}))

pdf(file = "./plots/figure_s6c.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = agglomerate_pc,
     pch = 16,
     col = "#00000040",
     cex = 0.75,
     xlim = c(-60,60),
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = NA,
     ylab = "Proliferation",
     las = 1,
     lwd = 0.5/0.75,axes = F)
lines(x = c(0,0),y = c(0,140),lwd = 0.5/0.75)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-10.5,2),lwd = 0.5/0.75)
lines(x = c(5,52),y = rep(-10.5,2),lwd = 0.5/0.75)
text(x = -27.5,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()

#reduced x range
pdf(file = "./plots/figure_s6c.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(1.6,0.6,0),xpd = T)
plot(x = agglomerate_pc,
     pch = 16,
     col = "#00000040",
     cex = 0.75,
     xlim = c(-60,60)/2,
     ylim = c(0,140),
     frame = F,
     xaxs = "i",
     yaxs = "i",
     xlab = NA,
     ylab = "Proliferation",
     las = 1,
     lwd = 0.5/0.75,axes = F)
lines(x = c(0,0),y = c(0,140),lwd = 0.5/0.75)
axis(side = 1,at = seq(-60,60,20)/2,labels = c(-60,NA,NA,0,NA,NA,60)/2,lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60/2,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143/2,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71)/2,y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52)/2,y = rep(-10.5,2),lwd = 0.5/0.75)
lines(x = c(5,52)/2,y = rep(-10.5,2),lwd = 0.5/0.75)
text(x = -27.5/2,y = -15,labels = "Proneural",cex = 6/7)
text(x = 27.5/2,y = -15,labels = "Mesenchymal",cex = 6/7)
dev.off()