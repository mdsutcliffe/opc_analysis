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


# Re-run to remove NFOs for mixture matrix
signature_cibersort <- signatureMatrix(geneList = bulk$rsem$symbol,removeNFO = T)

wt12_cibersort <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 12 & bulk$info$genotype == "WT"))]
wt12_cibersort <- wt12_cibersort[match(signature_cibersort$genes,wt12_cibersort$symbol),c(3,10:ncol(wt12_cibersort))]

wt150_cibersort <- bulk$tpm[,c(1:9,9+which(bulk$info$day == 150 & bulk$info$genotype == "WT"))]
wt150_cibersort <- wt150_cibersort[match(signature_cibersort$genes,wt150_cibersort$symbol),c(3,10:ncol(wt150_cibersort))]

write.table(x = wt12_cibersort,file = "./temp/mixture_bulk_12_wt.txt",quote = F,sep = "\t",row.names = F)
write.table(x = wt150_cibersort,file = "./temp/mixture_bulk_150_wt.txt",quote = F,sep = "\t",row.names = F)


# Plot
{
  pdf(file = "./plots/bulk_LCM_purity_150_WT.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mfrow = c(2,1),mgp = c(2.9,1,0))
  boxplot(t(bulk_tpm_150_WT_sig),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = "Bulk TPM / (Signature TPM)")
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  boxplot(t(bulk_tpm_150_WT_sig),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylim = c(0,0.3),ylab = "Bulk TPM / (Signature TPM)")
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  dev.off()
  
  pdf(file = "./plots/bulk_LCM_purity_12_WT.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mfrow = c(2,1),mgp = c(2.9,1,0))
  boxplot(t(bulk_tpm_12_WT_sig),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = "Bulk TPM / (Signature TPM)")
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  boxplot(t(bulk_tpm_12_WT_sig),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylim = c(0,0.3),ylab = "Bulk TPM / (Signature TPM)")
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  dev.off()
  
  pdf(file = "./plots/bulk_LCM_purity_150_WT_log2.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mgp = c(2.9,1,0))
  boxplot(t(bulk_tpm_150_WT_sig_log2),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = expression("Log"[2]*"(Bulk TPM) - Log"[2]*"(Signature TPM)"),ylim = c(-12,2))
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  dev.off()
  
  pdf(file = "./plots/bulk_LCM_purity_12_WT_log2.pdf",width = 10,height = 4)
  par(mar = c(6,4,1,1),mgp = c(2.9,1,0))
  boxplot(t(bulk_tpm_12_WT_sig_log2),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = expression("Log"[2]*"(Bulk TPM) - Log"[2]*"(Signature TPM)"),ylim = c(-12,2))
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  dev.off()
}





# What are their highest TPM values, what are our highest?
# bulk_tpm_150_WT_sort <- sort(x = rowSums(x = bulk_tpm_150_WT[,10:ncol(bulk_tpm_150_WT)]) / length(10:ncol(bulk_tpm_150_WT)),decreasing = T)

bulk_median <- apply(X = bulk_tpm_150_WT[,10:ncol(bulk_tpm_150_WT)],MARGIN = 1,FUN = median)
bulk_tpm_150_WT_sort <- bulk_tpm_150_WT[order(as.numeric(bulk_median),decreasing = T),]
bulk_tpm_150_WT_sort$opcBulk_median <- sort(bulk_median,decreasing = T)

sig_max <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 1,FUN = max)
sig_avg_tpm_sort <- sig_avg_tpm[order(as.numeric(sig_max),decreasing = T),]

bulk_top <- cbind(bulk_tpm_150_WT_sort[1:10,10:ncol(bulk_tpm_150_WT_sort)],sig_avg_tpm[match(bulk_tpm_150_WT_sort$symbol[1:10],sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)])
names(bulk_top)[1:length(10:ncol(bulk_tpm_150_WT))] <- paste("bulk_WT_150",1:length(10:ncol(bulk_tpm_150_WT)),sep = "_")
write.csv(x = bulk_top,file = "./temp/LCM_purity_mostAbundant_bulk.csv",quote = F,row.names = T)

sig_top <- cbind(bulk_tpm_150_WT_sort[match(sig_avg_tpm_sort$symbol[1:10],bulk_tpm_150_WT_sort$symbol),10:ncol(bulk_tpm_150_WT_sort)],sig_avg_tpm_sort[1:10,2:ncol(sig_avg_tpm_sort)])

# Collect all in table
x_summarize <- rbind(cbind(bulk_tpm_150_WT[match(opc_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],sig_avg_tpm[match(opc_markers,sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)]),
                     cbind(bulk_tpm_150_WT[match(nfo_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],sig_avg_tpm[match(nfo_markers,sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)]),
                     cbind(bulk_tpm_150_WT[match(mo_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],sig_avg_tpm[match(mo_markers,sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)]),
                     cbind(bulk_tpm_150_WT[match(microglia_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],sig_avg_tpm[match(microglia_markers,sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)]),
                     cbind(bulk_tpm_150_WT[match(neuron_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],sig_avg_tpm[match(neuron_markers,sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)]),
                     cbind(bulk_tpm_150_WT[match(astrocyte_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],sig_avg_tpm[match(astrocyte_markers,sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)]),
                     cbind(bulk_tpm_150_WT[match(endothelial_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],sig_avg_tpm[match(endothelial_markers,sig_avg_tpm$symbol),2:ncol(sig_avg_tpm)]))
names(x_summarize)[1:length(10:ncol(bulk_tpm_150_WT))] <- paste("bulk_WT_150",1:length(10:ncol(bulk_tpm_150_WT)),sep = "_")
write.csv(x = x_summarize,file = "./temp/LCM_purity.csv",quote = F,row.names = T)
