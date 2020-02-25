# Verify LCM purity
setwd("/Users/mdsutcliffe/Github/opc_analysis")

library(readxl)
library(pheatmap)
library(RColorBrewer)
library(grid)

source("./functions/normalizeTPM.R")

# Pre-process bulk
{
  # File paths
  bulk.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/rsem_opcBulk.csv"
  bulk.info.path <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/info_opcBulk.csv"
  
  # Import
  bulk <- read.csv(bulk.path,stringsAsFactors = F)
  bulk.info <- read.csv(bulk.info.path,stringsAsFactors = F)
  
  # Keep only tdT+ samples
  bulk <- bulk[,c(1:9,9+which(bulk.info$celltype == "tdTpositive"))]
  bulk.info <- bulk.info[bulk.info$celltype == "tdTpositive",]
  
  # Remove outliers
  outliers <- c("opcBulk_90dpi_CKO_20647_tdTpositive_GAGCT","opcBulk_90dpi_WT_21167_tdTpositive")
  bulk <- bulk[,!(names(bulk) %in% outliers)]
  bulk.info <- bulk.info[!(bulk.info$name %in% outliers),]
  
  # Rename samples if they're from the same mouse as an outlier
  for (iOutlier in outliers) {
    if (length(strsplit(x = iOutlier,split = "_")[[1]]) == 6) {
      iRootName <- paste(strsplit(x = iOutlier,split = "_")[[1]][-6],collapse = "_")
      names(bulk)[which(grepl(iRootName,names(bulk)))] <- iRootName
      bulk.info$name[which(grepl(iRootName,bulk.info$name))] <- iRootName
    }
  }
  
  bulk_collapse <- bulk
  bulk_collapse.info <- bulk.info
  
  # Convert all columns to characters
  bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],as.character)
  
  # Collapse replicates from same mouse
  bulk_collapse_unique_replicateIndex.info <- unique(bulk_collapse.info[,c("day","genotype","sex","brainID")])
  indsToRemove <- c()
  for (i in 1:nrow(bulk_collapse_unique_replicateIndex.info)) {
    inds <- which(bulk_collapse.info$day == bulk_collapse_unique_replicateIndex.info$day[i] &
                    bulk_collapse.info$genotype == bulk_collapse_unique_replicateIndex.info$genotype[i] &
                    bulk_collapse.info$sex == bulk_collapse_unique_replicateIndex.info$sex[i] &
                    bulk_collapse.info$brainID == bulk_collapse_unique_replicateIndex.info$brainID[i])
    if (length(inds) > 1) {
      bulk_collapse[,9+inds[1]] <- rowSums(bulk_collapse[,9+inds])
      bulk_collapse.info[inds[1],"nCells"] <- paste(bulk_collapse.info[inds,"nCells"],collapse="|")
      bulk_collapse.info[inds[1],"runID"] <- paste(bulk_collapse.info[inds,"runID"],collapse="|")
      bulk_collapse.info[inds[1],"barcode"] <- paste(bulk_collapse.info[inds,"barcode"],collapse="|")
      bulk_collapse.info[inds[1],"replicate"] <- paste(bulk_collapse.info[inds,"replicate"],collapse="|")
      indsToRemove <- c(indsToRemove,inds[2:length(inds)])
      
      iNewName <- paste(strsplit(x = names(bulk_collapse)[9+inds[1]],split = "_")[[1]][-6],collapse = "_")
      names(bulk_collapse)[9+inds[1]] <- iNewName
      bulk_collapse.info$name[inds[1]] <- iNewName
    }
  }
  bulk_collapse <- bulk_collapse[,-(9+indsToRemove)]
  bulk_collapse.info <- bulk_collapse.info[-indsToRemove,]
  
  # Convert all columns to factors
  bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],factor)
  
}

# Pre-process signature
{
  nGenesEach <- 5
  f.barres <- list.files(path = "./external/GSE52564_RAW",full.names = T)
  f.barres_base <- basename(f.barres)
  
  sig <- read_xls(path = f.barres[1])[,1]
  names(sig)[1] <- "symbol"
  sig_samples <- do.call(cbind,sapply(X = 1:length(f.barres),FUN = function(i) read_xls(path = f.barres[i])[,2]))
  
  # Get signature names
  f.barres_base_cut <- sapply(X = 1:length(f.barres_base),FUN = function(i) strsplit(x = f.barres_base[i],split = "[_.]")[[1]][2])
  sig <- cbind(sig,sig_samples)
  names(sig)[2:ncol(sig)] <- f.barres_base_cut
  
  # Remove whole cortex
  sig <- sig[,1:15]
  f.barres_base_cut <- f.barres_base_cut[1:14]
  
  # Remove NFO
  sig <- sig[,-c(7,8)]
  f.barres_base_cut <- f.barres_base_cut[-c(7,8)]
  
  # Average replicates
  cell_types <- sapply(X = f.barres_base_cut,FUN = function(x) substr(x = x,start = 1,stop = nchar(x)-1),USE.NAMES = F)
  cell_types <- factor(cell_types)
  sig_avg <- cbind(sig[,1,drop = F],do.call(cbind,lapply(X = 1:length(levels(cell_types)),FUN = function(i) rowMeans(sig[,1 + which(cell_types == levels(cell_types)[i])]))))
  names(sig_avg)[2:ncol(sig_avg)] <- levels(cell_types)
}

commonGenes <- intersect(x = sig_avg$symbol,y = bulk_collapse$symbol)

# Convert FPKM to TPM
sig_avg_tpm <- sig_avg[sig_avg$symbol %in% commonGenes,]
sig_avg_tpm[,2:ncol(sig_avg_tpm)] <- apply(X = sig_avg_tpm[,2:ncol(sig_avg_tpm)],MARGIN = 2,FUN = function(x) x / sum(x) * 10^6)

bulk_tpm <- normalizeTPM(rsem = bulk_collapse[bulk_collapse$symbol %in% commonGenes,],index_counts = 10:ncol(bulk_collapse))

# Get only 150 dpi WT
bulk_tpm_150_WT <- bulk_tpm[,c(1:9,9+which(bulk_collapse.info$day == 150 & bulk_collapse.info$genotype == "WT"))]
row.names(bulk_tpm_150_WT) <- bulk_tpm_150_WT$symbol

# Get only 12 dpi WT
bulk_tpm_12_WT <- bulk_tpm[,c(1:9,9+which(bulk_collapse.info$day == 12 & bulk_collapse.info$genotype == "WT"))]
row.names(bulk_tpm_12_WT) <- bulk_tpm_12_WT$symbol

opc_markers <- c("Pdgfra","Lnx1","Dcn","Mmp15","Cdo1")
nfo_markers <- c("Gp1bb","Tmem108","Fyn","Ust","Mical3")
mo_markers <- c("Gjb1","Ndrg1","Ppp1r14a","Adssl1","Aspa")
pericyte_markers <- c("Fmod","Rps2","Igf2","Gpc3","Ogn")
microglia_markers <- c("Slfn2","Gpr84","Ccr7","Bcl2a1d","Tnf")
neuron_markers <- c("Reln","Nhlh2","Slc17a6","Trp73","Lhx5")
astrocyte_markers <- c("Hgf","Aqp4","Itih3","Bmpr1b","Itga7")
endothelial_markers <- c("Cldn5","Ttr","Ly6a","Madcam1","Akr1c14")

# Divided
{
  bulk_tpm_150_WT_opc <- apply(X = bulk_tpm_150_WT[match(opc_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(opc_markers,sig_avg_tpm$symbol),"OPC"])
  
  bulk_tpm_150_WT_nfo <- apply(X = bulk_tpm_150_WT[match(nfo_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(nfo_markers,sig_avg_tpm$symbol),"NFO"])
  
  bulk_tpm_150_WT_mo <- apply(X = bulk_tpm_150_WT[match(mo_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(mo_markers,sig_avg_tpm$symbol),"MO"])
  
  bulk_tpm_150_WT_microglia <- apply(X = bulk_tpm_150_WT[match(microglia_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(microglia_markers,sig_avg_tpm$symbol),"Microglia"])
  
  bulk_tpm_150_WT_neuron <- apply(X = bulk_tpm_150_WT[match(neuron_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(neuron_markers,sig_avg_tpm$symbol),"Neuron"])
  
  bulk_tpm_150_WT_astrocyte <- apply(X = bulk_tpm_150_WT[match(astrocyte_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(astrocyte_markers,sig_avg_tpm$symbol),"Astrocyte"])
  
  bulk_tpm_150_WT_endothelial <- apply(X = bulk_tpm_150_WT[match(endothelial_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(endothelial_markers,sig_avg_tpm$symbol),"Endothelial"])
  
  bulk_tpm_150_WT_sig <- rbind(bulk_tpm_150_WT_opc,
                               bulk_tpm_150_WT_nfo,
                               bulk_tpm_150_WT_mo,
                               bulk_tpm_150_WT_microglia,
                               bulk_tpm_150_WT_neuron,
                               bulk_tpm_150_WT_astrocyte,
                               bulk_tpm_150_WT_endothelial)
  
  
  bulk_tpm_12_WT_opc <- apply(X = bulk_tpm_12_WT[match(opc_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(opc_markers,sig_avg_tpm$symbol),"OPC"])
  
  bulk_tpm_12_WT_nfo <- apply(X = bulk_tpm_12_WT[match(nfo_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(nfo_markers,sig_avg_tpm$symbol),"NFO"])
  
  bulk_tpm_12_WT_mo <- apply(X = bulk_tpm_12_WT[match(mo_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(mo_markers,sig_avg_tpm$symbol),"MO"])
  
  bulk_tpm_12_WT_microglia <- apply(X = bulk_tpm_12_WT[match(microglia_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(microglia_markers,sig_avg_tpm$symbol),"Microglia"])
  
  bulk_tpm_12_WT_neuron <- apply(X = bulk_tpm_12_WT[match(neuron_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(neuron_markers,sig_avg_tpm$symbol),"Neuron"])
  
  bulk_tpm_12_WT_astrocyte <- apply(X = bulk_tpm_12_WT[match(astrocyte_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(astrocyte_markers,sig_avg_tpm$symbol),"Astrocyte"])
  
  bulk_tpm_12_WT_endothelial <- apply(X = bulk_tpm_12_WT[match(endothelial_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) x / sig_avg_tpm[match(endothelial_markers,sig_avg_tpm$symbol),"Endothelial"])
  
  bulk_tpm_12_WT_sig <- rbind(bulk_tpm_12_WT_opc,
                              bulk_tpm_12_WT_nfo,
                              bulk_tpm_12_WT_mo,
                              bulk_tpm_12_WT_microglia,
                              bulk_tpm_12_WT_neuron,
                              bulk_tpm_12_WT_astrocyte,
                              bulk_tpm_12_WT_endothelial)
}

# Log2 transformed and subtracted
{
  bulk_tpm_150_WT_opc_log2 <- apply(X = bulk_tpm_150_WT[match(opc_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(opc_markers,sig_avg_tpm$symbol),"OPC"] + 1))
  
  bulk_tpm_150_WT_nfo_log2 <- apply(X = bulk_tpm_150_WT[match(nfo_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(nfo_markers,sig_avg_tpm$symbol),"NFO"] + 1))
  
  bulk_tpm_150_WT_mo_log2 <- apply(X = bulk_tpm_150_WT[match(mo_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(mo_markers,sig_avg_tpm$symbol),"MO"] + 1))
  
  bulk_tpm_150_WT_microglia_log2 <- apply(X = bulk_tpm_150_WT[match(microglia_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(microglia_markers,sig_avg_tpm$symbol),"Microglia"] + 1))
  
  bulk_tpm_150_WT_neuron_log2 <- apply(X = bulk_tpm_150_WT[match(neuron_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(neuron_markers,sig_avg_tpm$symbol),"Neuron"] + 1))
  
  bulk_tpm_150_WT_astrocyte_log2 <- apply(X = bulk_tpm_150_WT[match(astrocyte_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(astrocyte_markers,sig_avg_tpm$symbol),"Astrocyte"] + 1))
  
  bulk_tpm_150_WT_endothelial_log2 <- apply(X = bulk_tpm_150_WT[match(endothelial_markers,bulk_tpm_150_WT$symbol),10:ncol(bulk_tpm_150_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(endothelial_markers,sig_avg_tpm$symbol),"Endothelial"] + 1))
  
  bulk_tpm_150_WT_sig_log2 <- rbind(bulk_tpm_150_WT_opc_log2,
                                    bulk_tpm_150_WT_nfo_log2,
                                    bulk_tpm_150_WT_mo_log2,
                                    bulk_tpm_150_WT_microglia_log2,
                                    bulk_tpm_150_WT_neuron_log2,
                                    bulk_tpm_150_WT_astrocyte_log2,
                                    bulk_tpm_150_WT_endothelial_log2)
  
  bulk_tpm_12_WT_opc_log2 <- apply(X = bulk_tpm_12_WT[match(opc_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(opc_markers,sig_avg_tpm$symbol),"OPC"] + 1))
  
  bulk_tpm_12_WT_nfo_log2 <- apply(X = bulk_tpm_12_WT[match(nfo_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(nfo_markers,sig_avg_tpm$symbol),"NFO"] + 1))
  
  bulk_tpm_12_WT_mo_log2 <- apply(X = bulk_tpm_12_WT[match(mo_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(mo_markers,sig_avg_tpm$symbol),"MO"] + 1))
  
  bulk_tpm_12_WT_microglia_log2 <- apply(X = bulk_tpm_12_WT[match(microglia_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(microglia_markers,sig_avg_tpm$symbol),"Microglia"] + 1))
  
  bulk_tpm_12_WT_neuron_log2 <- apply(X = bulk_tpm_12_WT[match(neuron_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(neuron_markers,sig_avg_tpm$symbol),"Neuron"] + 1))
  
  bulk_tpm_12_WT_astrocyte_log2 <- apply(X = bulk_tpm_12_WT[match(astrocyte_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(astrocyte_markers,sig_avg_tpm$symbol),"Astrocyte"] + 1))
  
  bulk_tpm_12_WT_endothelial_log2 <- apply(X = bulk_tpm_12_WT[match(endothelial_markers,bulk_tpm_12_WT$symbol),10:ncol(bulk_tpm_12_WT)],MARGIN = 2,FUN = function(x) log2(x + 1) - log2(sig_avg_tpm[match(endothelial_markers,sig_avg_tpm$symbol),"Endothelial"] + 1))
  
  bulk_tpm_12_WT_sig_log2 <- rbind(bulk_tpm_12_WT_opc_log2,
                                   bulk_tpm_12_WT_nfo_log2,
                                   bulk_tpm_12_WT_mo_log2,
                                   bulk_tpm_12_WT_microglia_log2,
                                   bulk_tpm_12_WT_neuron_log2,
                                   bulk_tpm_12_WT_astrocyte_log2,
                                   bulk_tpm_12_WT_endothelial_log2)
}

# Cibersort
{
  sig_avg <- sig_avg[match(commonGenes,sig_avg$symbol),]
  
  sig_Genes <- lapply(X = levels(cell_types),FUN = function(i) {
    iSig <- sig_avg[,c(1,1+which(levels(cell_types) == i))]
    
    # if (i %in% c("OPC","NFO","MO")) {
    #   iOther <- sig_avg[,c(1,1+which(!(levels(cell_types) %in% c("OPC","NFO","MO"))))]
    # } else {
    iOther <- sig_avg[,c(1,1+which(levels(cell_types) != i))]
    # }
    
    iSig$fc <- iSig[,2] / rowMeans(iOther[,2:ncol(iOther)])
    
    iSig <- iSig[order(iSig$fc,decreasing = T),]
    iSig <- iSig[iSig[,2] > 20,]
    
    return(iSig$symbol[1:1000])
  })
  
  sig_mat <- do.call(cbind,sig_Genes)
  colnames(sig_mat) <- levels(cell_types)
  
  for (i in 2:nrow(sig_mat)) {
    for (j in 1:ncol(sig_mat)) {
      if (sig_mat[i,j] %in% sig_mat[1:(i-1),]) {
        sig_mat[i,j] = NA
      }
      if (j > 1 & sig_mat[i,j] %in% sig_mat[i,1:(j-1)]) {
        sig_mat[i,j] = NA
      }
    }
  }
  
  sig_Genes_unique <- apply(sig_mat,2,function(x) {
    y <- na.omit(x)
    return(y[1:40])
  })
  
  allSigGenes <- c()
  for (i in 1:ncol(sig_Genes_unique)) {
    allSigGenes <- c(allSigGenes,sig_Genes_unique[1:40,i])
  }
  
  bulk_tpm_cibersort_150_WT <- bulk_tpm_150_WT[match(allSigGenes,bulk_tpm_150_WT$symbol),c(3,10:ncol(bulk_tpm_150_WT))]
  write.table(x = bulk_tpm_cibersort_150_WT,file = "./temp/mixture_bulk150wt.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  
  sig_tpm_cibersort <- sig_avg_tpm[match(allSigGenes,sig_avg_tpm$symbol),c(2:ncol(sig_avg_tpm))]
  
  sig_mixture <- sig_avg_tpm
  names(sig_mixture)[2:ncol(sig_mixture)] <- paste0(names(sig_mixture)[2:ncol(sig_mixture)],"_mixture")
  write.table(x = sig_mixture,file = "./temp/mixture_signature.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  
  # results
  x <- read.table(file = "~/Downloads/CIBERSORTx_Job7_Results.txt",header = T,sep = "\t")
  row.names(x) <- x$Mixture
  y <- x[,2:7]
  x$Absolute.score..sig.score.
  pdf(file = "./plots/cibersort_bulk_150_wt_absolute.pdf",width = 8,height = 8)
  pheatmap(y,col = (brewer.pal(9,"Reds")),labels_row = rep("  ",nrow(y)))
  grid.text(label = "150 dpi WT bulk samples",x = 0.95,rot = 90)
  dev.off()
  
  z <- y
  z$Other <- 100 - rowSums(z)
  pdf(file = "./plots/cibersort_bulk_150_wt_other.pdf",width = 8,height = 8)
  pheatmap(z,col = (brewer.pal(9,"Reds")),labels_row = rep("  ",nrow(y)),breaks = seq(0,100,length.out = 10))
  grid.text(label = "150 dpi WT bulk samples",x = 0.95,rot = 90)
  dev.off()
  
  pdf(file = "./plots/cibersort_bulk_150_wt_absolute.pdf",width = 8,height = 8)
  pheatmap(z,col = (brewer.pal(9,"Reds")),labels_row = rep("  ",nrow(z)),breaks = seq(0,6,length.out = 10))
  grid.text(label = "150 dpi WT bulk samples",x = 0.95,rot = 90)
  dev.off()
  
  
  
  
  x <- read.table(file = "~/Downloads/CIBERSORTx_Job8_Results.txt",header = T,sep = "\t")
  row.names(x) <- x$Mixture
  y <- x[,2:7]
  pdf(file = "./plots/cibersort_signature_absolute.pdf",width = 5.8,height = 5)
  pheatmap(y,col = (brewer.pal(9,"Reds")),display_numbers = T,number_color = "#FFFFFF",breaks = seq(0,100,length.out = 10))
  dev.off()
}


# Plot
{
  pdf(file = "./plots/bulk_LCM_purity_150_WT.pdf",width = 10,height = 8)
  par(mar = c(6,4,1,1),mfrow = c(2,1),mgp = c(2.9,1,0))
  boxplot(t(bulk_tpm_150_WT_sig),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylab = "Bulk TPM / (Signature TPM)")
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  boxplot(t(bulk_tpm_150_WT_sig),at = (1:49)[(1:49) %% 7 != 0 & (1:49) %% 7 != 6],las = 2,ylim = c(0,0.3),ylab = "Bulk TPM / (Signature TPM)")
  title(xlab = "     OPC                    NFO                     MO                 Microglia                Neuron              Astrocyte               Endothelial",mgp = c(5,0,0))
  dev.off()
  
  pdf(file = "./plots/bulk_LCM_purity_12_WT.pdf",width = 10,height = 8)
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
