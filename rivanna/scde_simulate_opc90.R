library(scde)

args <- commandArgs(T)
iteration <- as.numeric(args[1])
message(paste("Iteration:",iteration))
set.seed(iteration)

f.opc90 <- "./rsem_opc90.csv"
f.opc90.info <- "./info_opc90.csv"

opc90 <- read.csv(file = f.opc90, stringsAsFactors = F)
opc90.info <- read.csv(file = f.opc90.info, stringsAsFactors = F)

iFemalePooled <- sample(which(opc90.info$sex == "female" & opc90.info$type == "pooled"))
iFemaleTenCell <- sample(which(opc90.info$sex == "female" & opc90.info$type == "ten-cell"))
iMalePooled <- sample(which(opc90.info$sex == "male" & opc90.info$type == "pooled"))
iMaleTenCell <- sample(which(opc90.info$sex == "male" & opc90.info$type == "ten-cell"))

dropAndDuplicate <- function(index) {
  index <- index[c(1:(length(index)-1),length(index)-1)]
}

iFemalePooled <- dropAndDuplicate(iFemalePooled)
iFemaleTenCell <- dropAndDuplicate(iFemaleTenCell)
iMalePooled <- dropAndDuplicate(iMalePooled)
iMaleTenCell <- dropAndDuplicate(iMaleTenCell)

runPagoda <- function(rsem,iPooled,iTenCell) {
  nPooled <- length(iPooled)
  nTenCell <- length(iTenCell)
  
  rsem <- rsem[rsem$chr != "ERCC" & rsem$chr != "MT",]
  row.names(rsem) <- rsem$symbol
  
  rsem <- rsem[,9+c(iPooled,iTenCell)]
  
  rsem_clean <- clean.counts(counts = rsem,min.lib.size = 1000,min.reads = 1,min.detected = 1)
  rsem_clean <- as.matrix(rsem_clean)
  mode(rsem_clean) <- "integer"
  
  rsem_clean_pooled <- rsem_clean[,1:length(iPooled)]
  rsem_clean_tencell <- rsem_clean[,(length(iPooled)+1):ncol(rsem_clean)]
  
  models_pooled <- knn.error.models(counts = rsem_clean_pooled,k = ncol(rsem_clean_pooled)/4,min.nonfailed = 5,min.count.threshold = 5,save.model.plots = F,n.cores = 1)
  models_tencell <- knn.error.models(counts = rsem_clean_tencell,k = ncol(rsem_clean_tencell)/4,min.nonfailed = 5,min.count.threshold = 5,save.model.plots = F,n.cores = 1)
  
  varinfo_pooled <- pagoda.varnorm(models = models_pooled,counts = rsem_clean_pooled,trim = 3/ncol(rsem_clean_pooled),plot = F,n.cores = 1)
  varinfo_tencell <- pagoda.varnorm(models = models_tencell,counts = rsem_clean_tencell,trim = 3/ncol(rsem_clean_tencell),plot = F,n.cores = 1)
  
  varinfo_pooled <- pagoda.subtract.aspect(varinfo = varinfo_pooled,aspect = colSums(rsem_clean_pooled[,row.names(models_pooled)] > 0))
  varinfo_tencell <- pagoda.subtract.aspect(varinfo = varinfo_tencell,aspect = colSums(rsem_clean_tencell[,row.names(models_tencell)] > 0))
  
  arv <- merge(x = data.frame(tencell_var = varinfo_tencell$arv),y = data.frame(pooled_var = varinfo_pooled$arv),by = "row.names")
  names(arv)[1] <- "symbol"
  
  threshold <- quantile(arv$pooled_var,0.95)
  
  candidates <- as.character(arv$symbol[arv$tencell_var > threshold & arv$pooled_var <= threshold])
  
  return(list(candidates = candidates,
              threshold = threshold,
              arv = arv,
              varinfo_tencell = varinfo_tencell,
              varinfo_pooled = varinfo_pooled))
}
message("\nStarting Female...")
message(paste("Pooled:  ",paste(sprintf("%02d",iFemalePooled),collapse = " ")))
message(paste("Samples: ",paste(sprintf("%02d",iFemaleTenCell),collapse = " ")))

resFemale <- runPagoda(rsem = opc90,
                       iPooled = iFemalePooled,
                       iTenCell = iFemaleTenCell)
message("Finished Female!")

message("\nStarting Male...")
message(paste("Pooled:  ",paste(sprintf("%02d",iMalePooled),collapse = " ")))
message(paste("Samples: ",paste(sprintf("%02d",iMaleTenCell),collapse = " ")))

resMale <- runPagoda(rsem = opc90,
                       iPooled = iMalePooled,
                       iTenCell = iMaleTenCell)
message("Finished Male!")

message("\nWriting files...")
path_RData <- "./rdatafiles"
path_candidates <- "./candidates"
path_arv <- "./arv"

if (!dir.exists(paths = path_RData)) dir.create(path = path_RData)
save.image(file = paste0(path_RData,"/simulation_",sprintf("%03d",iteration),".RData"))

if (!dir.exists(paths = path_candidates)) dir.create(path = path_candidates)
write.table(x = resFemale$candidates,file = paste0(path_candidates,"/candidates_female_",sprintf("%03d",iteration),".tsv"),quote = F,sep = "\t",row.names = F,col.names = F)
write.table(x = resMale$candidates,file = paste0(path_candidates,"/candidates_male_",sprintf("%03d",iteration),".tsv"),quote = F,sep = "\t",row.names = F,col.names = F)

if (!dir.exists(paths = path_arv)) dir.create(path = path_arv)
write.table(x = resFemale$arv,file = paste0(path_arv,"/arv_female_",sprintf("%03d",iteration),".tsv"),quote = F,sep = "\t",row.names = F)
write.table(x = resMale$arv,file = paste0(path_arv,"/arv_male_",sprintf("%03d",iteration),".tsv"),quote = F,sep = "\t",row.names = F)

message("\nCompleted successfully!")
