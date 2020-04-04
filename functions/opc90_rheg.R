
source("./import/import_opc90.R")

opc90_rheg <- function() {
  resF <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc90_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
  resM <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc90_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
  
  nAppearancesF <- table(unlist(x = resF))
  nAppearancesM <- table(unlist(x = resM))
  
  opc90_uniqueF <- setdiff(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])
  opc90_uniqueM <- setdiff(x = names(nAppearancesM)[nAppearancesM >= 75],y = names(nAppearancesF)[nAppearancesF >= 75])
  opc90_rheg <- intersect(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])
  
  return(list(Fspecific = opc90_uniqueF,
              Mspecific = opc90_uniqueM,
              RHEG = opc90_rheg))
}

opc90$candidates <- opc90_rheg()