library(scde)

source("./runPagoda.R")

load("./opc12_prepared.RData")

args <- commandArgs(T)
iteration <- as.numeric(args[1])
print(paste("Iteration:",iteration))

leaveOneOut <- function(rsem,iControls,iSamples,save_path) {
  
  print(paste("Controls:",paste(sprintf("%02d",iControls),collapse = " ")))
  print(paste("Samples: ",paste(sprintf("%02d",iSamples),collapse = " ")))
  
  res <- runPagoda(rsem = rsem,
                   iControls = iControls,
                   iSamples = iSamples)
  
  df <- rsem[which(rsem$symbol %in% res$genes),]
  
  if (!dir.exists(paths = save_path)) {
    dir.create(path = save_path)
  }
  write.table(x = df,
              quote = F,
              sep = "\t",
              row.names = F,
              file = paste0(save_path,"/resample_",sprintf("%03d",iteration),".tsv"))
  
  write.table(x = res$adjusted_variance,
              quote = F,
              sep = "\t",
              row.names = F,
              file = paste0(save_path,"/resample_",sprintf("%03d",iteration),"_variance.tsv"))
}

leaveOneOut(rsem = opc12_clean,
            iControls = iFemale_pooled[,iteration],
            iSamples = iFemale_tencell[,iteration],
            save_path <- "./results_female")

leaveOneOut(rsem = opc12_clean,
            iControls = iMale_pooled[,iteration],
            iSamples = iMale_tencell[,iteration],
            save_path <- "./results_male")
