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
plot(x = NA,y = NA,
     xlim = c(-60,60),ylim = c(0,8),
     xaxs = "i",yaxs = "i",
     xlab = NA,ylab = "Frequency",main = NA,
     frame = F,axes = F)
lines(x = c(0,0),y = c(0,7.5),lwd = 0.5/0.75)
hist(x = opc12_projection1[opc12$info$type == "ten-cell"],
     breaks = seq(-30,30,length.out = 23),
     col = "#969696",
     border = NA, add = T)
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
# Maximum likelihood estimate (MLE):
#              p_1             mu_1             mu_2            sigma 
#           0.8452           2.5930           1.5350           0.0790 
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
# p_1              0.80169854 0.8805800
# mu_1_gene_Gene 1 2.57343850 2.6125615
# mu_2_gene_Gene 1 1.34173518 1.7282648
# sigma            0.05483214 0.1138201
set.seed(0)
y <- replicate(n = 1e6,expr = {
  z <- replicate(n = 10,expr = {
    if (runif(n = 1) < ml$mle[1]) {
      return(rnorm(n = 1,mean = ml$mle[2],sd = ml$mle[4]))
    } else {
      return(rnorm(n = 1,mean = ml$mle[3],sd = ml$mle[4]))
    }
  })
  return(sum(exp(z)))
})

set.seed(0)
y_ci <- replicate(n = 1e6,expr = {
  z <- replicate(n = 10,expr = {
    if (runif(n = 1) < runif(n = 1,min = ml$ci[1,1],max = ml$ci[1,2])) {
      return(rnorm(n = 1,
                   mean = runif(n = 1,min = ml$ci[2,1],max = ml$ci[2,2]),
                   sd = runif(n = 1,min = ml$ci[4,1],max = ml$ci[4,2])))
    } else {
      return(rnorm(n = 1,
                   mean = runif(n = 1,min = ml$ci[3,1],max = ml$ci[3,2]),
                   sd = runif(n = 1,min = ml$ci[4,1],max = ml$ci[4,2])))
    }
  })
  return(sum(exp(z)))
})

save.image("./build/mle_results.RData")
load("./build/mle_results.RData")

step <- 0.5
x_range <- seq(-60,60,step)
ci <- sapply(X = x_range,FUN = function(i) {
  ivals <- (y_ci-shift)[(y_ci-shift) >= (i-step/2) & (y_ci-shift) < (i+step/2)]
  if (length(ivals) != 0) {
    return(diff(range(ivals))/2)
  } else {
    return(0)
  }
  
})

density_function <- approxfun(density(y-shift))
y2 <- density_function(x_range)
y2[is.na(y2)] <- 0

pdf(file = "./plots/figure_s6b.pdf",width = 2.25,height = 2.25,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0,0),mgp = c(2.6,0.6,0),xpd = T)
plot(x = NA,y = NA,
     xlim = c(-60,60),ylim = c(0,0.05),
     main = NA,xlab = NA,ylab = "Density",
     frame = F,axes = F,
     xaxs = "i",yaxs = "i",
     lwd = 0.5/0.75)
# polygon(x = c(-60:60,60:-60),y = c(y2 - ci*y2,rev(y2 + ci*y2)),border = NA,col = "#bdbdbd")
# polygon(x = c(seq(-60,60,0.1),rev(seq(-60,60,0.1))),y = c(y2_ci - ci_0.5*y2_ci,rev(y2_ci + ci_0.5*y2_ci)),border = NA,col = "#bdbdbd")
# polygon(x = c(x_range,rev(x_range)),y = c(y2-ci*y2,rev(y2+ci*y2)))
lines(x = c(0,0),y = c(0,(15/16)*0.05),lwd = 0.5/0.75)
df <- approxfun(density(y-shift))
df_x <- seq(-60,60,0.01)
df_y <- df(df_x)
df_y[is.na(df_y)] <- 0
lines(x = df_x,y = df_y,lwd = 1/0.75)
axis(side = 1,at = seq(-60,60,20),labels = c(-60,NA,NA,0,NA,NA,60),lwd = 0.5/0.75)
axis(side = 2,lwd = 0.5/0.75,las = 1)
lines(x = c(-5,-52),y = rep(-10.5,2)/140*0.05,lwd = 0.5/0.75)
lines(x = c(5,52),y = rep(-10.5,2)/140*0.05,lwd = 0.5/0.75)
text(x = -27.5,y = -15/140*0.05,labels = "Proneural",cex = 6/7)
text(x = 27.5,y = -15/140*0.05,labels = "Mesenchymal",cex = 6/7)
dev.off()

set.seed(0)
y_lower <- replicate(n = 1e6,expr = {
  z <- replicate(n = 10,expr = {
    if (runif(n = 1) < ml$ci[1,1]) {
      return(rnorm(n = 1,
                   mean = ml$ci[2,1],
                   sd = ml$ci[4,1]))
    } else {
      return(rnorm(n = 1,
                   mean = ml$ci[3,1],
                   sd = ml$ci[4,1]))
    }
  })
  return(sum(exp(z)))
})
set.seed(0)
y_upper <- replicate(n = 1e6,expr = {
  z <- replicate(n = 10,expr = {
    if (runif(n = 1) < ml$ci[1,2]) {
      return(rnorm(n = 1,
                   mean = ml$ci[2,2],
                   sd = ml$ci[4,2]))
    } else {
      return(rnorm(n = 1,
                   mean = ml$ci[3,2],
                   sd = ml$ci[4,2]))
    }
  })
  return(sum(exp(z)))
})
save.image("./build/mle_ci.RData")


pdf(file = "./plots/mle_ci.pdf",width = 4,height = 4)
par(mai = c(0.5,0.5,0,0))
plot(x = NA,y = NA,
     xlim = c(-60,60),ylim = c(0,0.06),
     main = NA,xlab = NA,ylab = "Density",
     frame = F,
     xaxs = "i",yaxs = "i",
     lwd = 0.5/0.75)
lines(density(y_lower-shift),col = "#e41a1c",lwd = 1)
lines(density(y_upper-shift),col = "#377eb8",lwd = 1)
lines(density(y-shift),lwd = 2)
legend(x = "topright",legend = c("all_lower","all_upper","mle"),lty = 1,col = c("#e41a1c","#377eb8","#000000"),lwd = 2)
dev.off()
#CI IS VERY LARGE
# START HERE

x_range <- seq(0,4,0.001)
z <- function(j) {
  x <- sapply(X = 1:10,FUN = function(c) {
    if (c <= j) {
      dnorm(x = x_range,mean = ml$mle[2],sd = ml$mle[4])
    } else {
      dnorm(x = x_range,mean = ml$mle[3],sd = ml$mle[4])
    }
  })
}

zz <- rowSums(sapply(X = 0:10,FUN = function(j) {
  choose(10,j) * ml$mle[1]^j * (1-ml$mle[1])^(10-j) * rowSums(z(j))
}))


pdf(file = "./plots/figure_s6b_inset.pdf",width = 0.6,height = 0.6,family = "ArialMT",pointsize = 7,useDingbats = F)
par(mai = c(0,0,0,0))
plot(x = NA,y = NA,
     xlim = c(-90,90),ylim = c(0,50),
     xlab = "",ylab = "",
     # xaxs = "i",yaxs = "i",
     lwd = 0.5/0.75,axes = F)
lines(x = c(0,0),y = c(0,(15/16)*50),lwd = 0.5/0.75)
lines(exp(x_range)*10-shift,zz,lwd = 1/0.75,xpd = T)
box(which = "plot",lwd = 0.5/0.75)
dev.off()




z_lower <- function(j) {
  x <- sapply(X = 1:10,FUN = function(c) {
    if (c <= j) {
      dnorm(x = x_range,mean = ml$ci[2,1],sd = ml$ci[4,1])
    } else {
      dnorm(x = x_range,mean = ml$ci[3,1],sd = ml$ci[4,1])
    }
  })
}

zz_lower <- rowSums(sapply(X = 0:10,FUN = function(j) {
  choose(10,j) * ml$ci[1,1]^j * (1-ml$ci[1,1])^(10-j) * rowSums(z_lower(j))
}))

z_upper <- function(j) {
  x <- sapply(X = 1:10,FUN = function(c) {
    if (c <= j) {
      dnorm(x = x_range,mean = ml$ci[2,2],sd = ml$ci[4,2])
    } else {
      dnorm(x = x_range,mean = ml$ci[3,2],sd = ml$ci[4,2])
    }
  })
}

zz_upper <- rowSums(sapply(X = 0:10,FUN = function(j) {
  choose(10,j) * ml$ci[1,2]^j * (1-ml$ci[1,2])^(10-j) * rowSums(z_upper(j))
}))

plot(exp(x_range)*10-shift,zz,lwd = 0.5/0.75,xpd = T,type = "l")
lines(exp(x_range)*10-shift,zz_lower,col = "red")
lines(exp(x_range)*10-shift,zz_upper,col = "red")

plot(x = NA,y = NA,
     xlim = c(-100,100),ylim = c(0,50))
lines(exp(x_range)*10-shift,zz_lower,col = "red")
lines(exp(x_range)*10-shift,zz_upper,col = "blue")
lines(exp(x_range)*10-shift,zz,lwd = 2)






sapply(X = 1:10,)

m_ <- sapply(X = seq(-60,60,0.1),FUN = function(i) {
  ivals <- (y_ci-shift)[(y_ci-shift) >= i-0.05 & (y_ci-shift) < (i+0.05)]
  if (length(ivals) != 0) {
    return(median(ivals))
  } else {
    return(0)
  }
  
})
y2_ <- df(seq(-60,60,0.1)); y2_[is.na(y2_)] <- 0
lines(seq(-60,60,0.1),m_/max(m_)*0.05)

lines(density(m_))






















ci_0.5 <- sapply(X = seq(-60,60,0.1),FUN = function(i) {
  ivals <- (y_ci-shift)[(y_ci-shift) >= i-0.05 & (y_ci-shift) < (i+0.05)]
  if (length(ivals) != 0) {
    return(diff(range(ivals))/2)
  } else {
    return(0)
  }
  
})
y2_ci <- df(seq(-60,60,0.1)); y2_ci[is.na(y2_ci)] <- 0

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

xvals <- seq(-60,60,0.1)
y1 <- dnorm(x = xvals,mean = ml$mle[3] - log(shift),sd = exp(ml$mle[4])) * (1 - ml$mle[1])
y2 <- dnorm(x = xvals,mean = ml$mle[2] - log(shift),sd = exp(ml$mle[4])) * ml$mle[1]
plot(x = NA,y = NA,
     xlim = range(xvals),ylim = c(0,1))
lines(x = xvals,y = y1)
lines(x = xvals,y = y2)

library(mixtools)
pc1.k2 <- normalmixEM(x = opc12_projection1+shift,k = 2)
pc1.k2$mu - shift

xvals <- seq(-60,60,0.1)
y1 <- dnorm(x = xvals,mean = pc1.k2$mu[1] - shift,sd = log(pc1.k2$sigma[1])) * (1 - ml$mle[1])
y2 <- dnorm(x = xvals,mean = pc1.k2$mu[2] - shift,sd = log(pc1.k2$sigma[2])) * ml$mle[1]
plot(x = NA,y = NA,
     xlim = range(xvals),ylim = c(0,.2))
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





stochprof.loop(model = "LN-LN",
               dataset = as.matrix(x),
               n = 10,
               TY = 2,
               genenames = NULL,
               par.range = rbind(c(0.2,0.8),c(shift,shift+60),c(shift-60,shift),c(0,60)),
               loops = 20)
