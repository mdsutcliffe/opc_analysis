
varF <- read.table(file = "./data/opc12F_variance.tsv",header = T,stringsAsFactors = F)
varM <- read.table(file = "./data/opc12M_variance.tsv",header = T,stringsAsFactors = F)

varF_candidate <- varF[varF$candidate,1]
varM_candidate <- varM[varM$candidate,1]

varF_candidate_table <- table(varF_candidate)
varM_candidate_table <- table(varM_candidate)

candidates_F <- names(varF_candidate_table[varF_candidate_table >= 75])
candidates_M <- names(varM_candidate_table[varM_candidate_table >= 75])

newcan <- intersect(candidates_F,candidates_M)


length(intersect(hetGenes_opc12_intersect,newcan))
