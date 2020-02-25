# Verify LCM purity by looking at the top 5 genes

library(readxl)

source("./opcBulk/import_opcBulk.R")
source("./functions/signatureMatrix.R")
source("./functions/normalizeTPM.R")

signature <- signatureMatrix(geneList = bulk$rsem$symbol,removeNFO = F)
commonGenes <- intersect(bulk$tpm$symbol,signature$tpm$symbol)

bulk$tpm <- normalizeTPM(rsem = bulk$rsem[match(commonGenes,bulk$rsem$symbol),],index_counts = 10:ncol(bulk$rsem))

wt12 <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 12 & bulk$info$genotype == "WT"))]
wt150 <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 150 & bulk$info$genotype == "WT"))]

marker_opc <- c("Pdgfra","Lnx1","Dcn","Mmp15","Cdo1")
marker_nfo <- c("Gp1bb","Tmem108","Fyn","Ust","Mical3")
marker_mo <- c("Gjb1","Ndrg1","Ppp1r14a","Adssl1","Aspa")
marker_microglia <- c("Slfn2","Gpr84","Ccr7","Bcl2a1d","Tnf")
marker_neuron <- c("Reln","Nhlh2","Slc17a6","Trp73","Lhx5")
marker_astrocyte <- c("Hgf","Aqp4","Itih3","Bmpr1b","Itga7")
marker_endothelial <- c("Cldn5","Ttr","Ly6a","Madcam1","Akr1c14")
# marker_pericyte <- c("Fmod","Rps2","Igf2","Gpc3","Ogn")

markerOrder <- c("OPC","NFO","MO","Microglia","Neuron","Astrocyte","Endothelial")
markers <- c(marker_opc,marker_nfo,marker_mo,marker_microglia,marker_neuron,marker_astrocyte,marker_endothelial)

wt12 <- wt12[match(markers,wt12$symbol),]
wt150 <- wt150[match(markers,wt150$symbol),]

sig <- signature$tpm[match(markers,signature$tpm$symbol),c("symbol",markerOrder)]
sig$value <- sapply(X = 1:length(markers),FUN = function(x) sig[x,floor((x - 1) / 5) + 1 + 1])
sig <- sig[,c("symbol","value")]

tpm12 <- cbind(data.frame(row.names = wt12$symbol),wt12[,10:ncol(wt12)])
tpm150 <- cbind(data.frame(row.names = wt150$symbol),wt150[,10:ncol(wt150)])

ratio12 <- cbind(data.frame(row.names = wt12$symbol),wt12[,10:ncol(wt12)] / sig$value)
ratio150 <- cbind(data.frame(row.names = wt150$symbol),wt150[,10:ncol(wt150)] / sig$value)

log2diff12 <- cbind(data.frame(row.names = wt12$symbol),log2(wt12[,10:ncol(wt12)] + 1) - log2(sig$value + 1))
log2diff150 <- cbind(data.frame(row.names = wt150$symbol),log2(wt150[,10:ncol(wt150)] + 1) - log2(sig$value + 1))

# TPM boxplots
{
  pdf(file = "./plots/LCM_purity_tpm_12_WT_boxplot.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mgp = c(2.9,1,0))
  boxplot(t(tpm12),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = "Bulk TPM",ylim = c(0,400))
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  dev.off()
  
  pdf(file = "./plots/LCM_purity_tpm_150_WT_boxplot.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mgp = c(2.9,1,0))
  boxplot(t(tpm150),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = "Bulk TPM",ylim = c(0,700))
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  dev.off()
}

# Ratio boxplots
{
  pdf(file = "./plots/LCM_purity_ratio_12_WT_boxplot.pdf",width = 10,height = 8)
  par(mar = c(6,4,1,1),mfrow = c(2,1),mgp = c(2.9,1,0))
  boxplot(t(ratio12),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = "Bulk TPM / Signature TPM")
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  boxplot(t(ratio12),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylim = c(0,0.3),ylab = "Bulk TPM / Signature TPM")
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  dev.off()
  
  pdf(file = "./plots/LCM_purity_ratio_150_WT_boxplot.pdf",width = 10,height = 8)
  par(mar = c(6,4,1,1),mfrow = c(2,1),mgp = c(2.9,1,0))
  boxplot(t(ratio150),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = "Bulk TPM / Signature TPM")
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  boxplot(t(ratio150),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylim = c(0,0.3),ylab = "Bulk TPM / Signature TPM")
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  dev.off()
}

# Log2diff boxplots
{
  pdf(file = "./plots/LCM_purity_log2diff_12_WT_boxplot.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mgp = c(2.9,1,0))
  boxplot(t(log2diff12),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = expression("Log"[2]*"(Bulk TPM) - Log"[2]*"(Signature TPM)"))
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  dev.off()
  
  pdf(file = "./plots/LCM_purity_log2diff_150_WT_boxplot.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mgp = c(2.9,1,0))
  boxplot(t(ratio150),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = expression("Log"[2]*"(Bulk TPM) - Log"[2]*"(Signature TPM)"))
  mtext(text = markerOrder,side = 1,at = seq(3,49,7),line = 5)
  dev.off()
}
