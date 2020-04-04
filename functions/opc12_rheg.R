
source("./import/import_opc12.R")

opc12_rheg <- function() {
  resF <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_female_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
  resM <- lapply(X = 0:99,FUN = function(x) scan(file = paste0("./data/opc12_candidates/candidates_male_",sprintf("%03d",x),".tsv"),what = "character",quiet = T))
  
  nAppearancesF <- table(unlist(x = resF))
  nAppearancesM <- table(unlist(x = resM))
  
  opc12_uniqueF <- setdiff(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])
  opc12_uniqueM <- setdiff(x = names(nAppearancesM)[nAppearancesM >= 75],y = names(nAppearancesF)[nAppearancesF >= 75])
  opc12_rheg <- intersect(x = names(nAppearancesF)[nAppearancesF >= 75],y = names(nAppearancesM)[nAppearancesM >= 75])
  
  return(list(Fspecific = opc12_uniqueF,
              Mspecific = opc12_uniqueM,
              RHEG = opc12_rheg))
}

opc12$candidates <- opc12_rheg()