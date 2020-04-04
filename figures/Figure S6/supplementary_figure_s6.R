# Supplementary Figure S6

library(numDeriv)
library(MASS)

source("./import/import_opc12.R")
source("./functions/pca_subtype.R")
f.stochprof <- list.files("./external/Stochprof/R",full.names = T)
for (i in f.stochprof) { source(i) }

opc12$pc <- pca_subtype(opc12)

pdf(file = "./Figure S6/figure_s6a.pdf",width = 3.125,height = 3.125,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = NA,y = NA,
     xlim = c(-60,60),ylim = c(0,8),
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Frequency",main = "10cRNA-seq 12 dpi samples",font.main = 1,cex.main = 8/7,
     axes = F)
lines(x = c(0,0),y = c(0,(15/16)*8),lwd = 0.5/0.75)
hist(x = opc12$pc$projection1[opc12$info$type == "ten-cell"],
     breaks = seq(-30,30,length.out = 23),
     col = "#969696",
     border = NA,
     add = T)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75,mgp = c(1.6,0.6,0))
axis(side = 2,at = 0:8,lwd = 0.5/0.75,mgp = c(1.6,0.6,0),las = 1)
lines(x = c(-5.1,-52),y = rep(-0.5,2),lwd = 0.5/0.75)
lines(x = c(5.1,52),y = rep(-0.5,2),lwd = 0.5/0.75)
text(x = -27.5,y = -0.75,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -0.75,labels = "Mesenchymal",cex = 6/7)
dev.off()

shift <- 120
opc12_pc1_shift <- as.numeric(opc12$pc$projection1[opc12$info$type == "ten-cell"]) + shift
ml <- stochasticProfilingML()
# 3 opc12_pc1_shift 1 no

# Maximum likelihood estimate (MLE):
#              p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
#           0.8453           2.5930           1.5350           0.0790 
# 
# Value of negative log-likelihood function at MLE:
# 212.2336 
# 
# Violation of constraints:
# none
# 
# BIC:
# 440.5686 
# 
# Approx. 95% confidence intervals for MLE:
#                       lower     upper
# p_1              0.80180149 0.8806727
# mu_1_gene_Gene 1 2.57343689 2.6125631
# mu_2_gene_Gene 1 1.34178205 1.7282179
# sigma            0.05483139 0.1138217

set.seed(0)
ml_pc1 <- replicate(n = 1e6,expr = {
  z <- replicate(n = 10,expr = {
    if (runif(n = 1) < ml$mle[1]) {
      return(rnorm(n = 1,mean = ml$mle[2],sd = ml$mle[4]))
    } else {
      return(rnorm(n = 1,mean = ml$mle[3],sd = ml$mle[4]))
    }
  })
  return(sum(exp(z)))
})

density_func <- approxfun(density(ml_pc1-shift))
density_x <- seq(-60,60,0.01)
density_y <- density_func(density_x)
density_y[is.na(density_y)] <- 0

pdf(file = "./Figure S6/figure_s6b.pdf",width = 3.125,height = 3.125,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(2.35,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = NA,y = NA,
     xlim = c(-60,60),ylim = c(0,0.05),
     xlab = NA,ylab = "Density",main = "Maximum-likelihood inference",font.main = 1,cex.main = 8/7,
     frame = F,axes = F,
     xaxs = "i",yaxs = "i",
     lwd = 0.5/0.75)
lines(x = c(0,0),y = c(0,(15/16)*0.05),lwd = 0.5/0.75)
lines(x = density_x,y = density_y,lwd = 1/0.75)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,lwd = 0.5/0.75,las = 1)
lines(x = c(-5.1,-52),y = rep(-0.5,2)/8*0.05,lwd = 0.5/0.75)
lines(x = c(5.1,52),y = rep(-0.5,2)/8*0.05,lwd = 0.5/0.75)
text(x = -27.5,y = -0.75/8*0.05,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -0.75/8*0.05,labels = "Mesenchymal",cex = 6/7)
dev.off()

x_ml <- seq(0,4,0.001)
z <- function(j) {
  x <- sapply(X = 1:10,FUN = function(c) {
    if (c <= j) {
      dnorm(x = x_ml,mean = ml$mle[2],sd = ml$mle[4])
    } else {
      dnorm(x = x_ml,mean = ml$mle[3],sd = ml$mle[4])
    }
  })
}

y_ml <- rowSums(sapply(X = 0:10,FUN = function(j) {
  choose(10,j) * ml$mle[1]^j * (1-ml$mle[1])^(10-j) * rowSums(z(j))
}))

pdf(file = "./Figure S6/figure_s6b_inset.pdf",width = 0.6,height = 0.6,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0,0,0,0))
plot(x = NA,y = NA,
     xlim = c(-90,90),ylim = c(0,50),
     xlab = "",ylab = "",
     lwd = 0.5/0.75,axes = F)
lines(x = c(0,0),y = c(0,(15/16)*50),lwd = 0.5/0.75)
lines(exp(x_ml)*10-shift,y_ml,lwd = 1/0.75,xpd = T)
box(which = "plot",lwd = 0.5/0.75)
dev.off()


# Agglomerations
set.seed(0)
nCells <- sample(x = seq(700,1000,10),size = 100,replace = T)
agglomerate_pc <- sapply(X = nCells,FUN = function(nC) {
  cellGroups <- sample(x = 9 + which(opc12$info$type == "ten-cell"),size = nC / 10,replace = T)
  
  iOPC12 <- data.frame(agglom = rowMeans(opc12$tpm[,cellGroups]))
  
  return(iOPC12)
})
opc12_agglomerate <- list(tpm = cbind(opc12$tpm[,1:9],do.call(what = cbind,args = agglomerate_pc)))
opc12_agglomerate$pc <- pca_subtype(opc12_agglomerate)

pdf(file = "./Figure S6/figure_s6c.pdf",width = 3.125,height = 3.125,bg = "white",pointsize = 7,useDingbats = F)
par(mai = c(0.75,0.5,0.25,0.5),mgp = c(2.35,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = opc12_agglomerate$pc$projection1,
     y = opc12_agglomerate$pc$projection2,
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
lines(x = c(0,0),y = c(0,140)*15/16,lwd = 0.5/0.75)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,at = seq(0,140,20),labels = c(0,rep(NA,6),140),lwd = 0.5/0.75,las = 1)
text(x = 60,y = 2.25,labels = "PC1",adj = c(0.5,0))
text(x = -58.07143,y = 140,labels = "PC2",adj = c(0,0.5))
polygon(x = c(-66,-66,-71),y = c(10,130,130),col = "#000000",border = NA)
lines(x = c(-5,-52),y = rep(-8.5,2),lwd = 0.5/0.75)
lines(x = c(5,52),y = rep(-8.5,2),lwd = 0.5/0.75)
text(x = -27.5,y = -12.75,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -12.75,labels = "Mesenchymal",cex = 6/7)
dev.off()