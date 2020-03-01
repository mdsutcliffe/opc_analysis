# PMID 30602536 - sex-specific gene cluster analysis

source("./opc12/RHEGs_opc12.R")

load("./build/orthologs.RData")

pmid30602536_table_S2 <- read.csv(file = "./external/PMID30602536_sex_specific_genes.csv",row.names = 1)

cluster_F_specific <- row.names(pmid30602536_table_S2)[pmid30602536_table_S2$female_cluster_genes & !pmid30602536_table_S2$male_cluster_genes]
cluster_M_specific <- row.names(pmid30602536_table_S2)[pmid30602536_table_S2$male_cluster_genes & !pmid30602536_table_S2$female_cluster_genes]
cluster_intersect <- row.names(pmid30602536_table_S2)[pmid30602536_table_S2$male_cluster_genes & pmid30602536_table_S2$female_cluster_genes]

opc12_uniqueF_human <- orthologs$HGNC.symbol[match(x = opc12_uniqueF,table = orthologs$MGI.symbol)]
opc12_uniqueM_human <- orthologs$HGNC.symbol[match(x = opc12_uniqueM,table = orthologs$MGI.symbol)]

opc12_uniqueF_human <- opc12_uniqueF_human[!is.na(opc12_uniqueF_human)]
opc12_uniqueM_human <- opc12_uniqueM_human[!is.na(opc12_uniqueM_human)]

orthologs$HGNC.symbol[!duplicated(orthologs$MGI.symbol) & !duplicated(orthologs$HGNC.symbol)]

F_intersect <- intersect(opc12_uniqueF_human,cluster_F_specific)
F_opc12_only <- setdiff(opc12_uniqueF_human,cluster_F_specific)
F_cluster_only <- setdiff(cluster_F_specific,opc12_uniqueF_human)
F_none <- setdiff(unique(orthologs$HGNC.symbol),union(opc12_uniqueF_human,cluster_F_specific))

fisher.test(x = matrix(data = c(length(F_intersect),length(F_opc12_only),
                                length(F_cluster_only),length(F_none)),
                       nrow = 2,ncol = 2))


M_intersect <- intersect(opc12_uniqueM_human,cluster_M_specific)
M_opc12_only <- setdiff(opc12_uniqueM_human,cluster_M_specific)
M_cluster_only <- setdiff(cluster_M_specific,opc12_uniqueM_human)
M_none <- setdiff(unique(orthologs$HGNC.symbol),union(opc12_uniqueM_human,cluster_M_specific))

fisher.test(x = matrix(data = c(length(M_intersect),length(M_opc12_only),
                                length(M_cluster_only),length(M_none)),
                       nrow = 2,ncol = 2))

cluster_F_specific_mouse <- orthologs$MGI.symbol[match(x = cluster_F_specific,table = orthologs$HGNC.symbol)]
cluster_F_specific_mouse <- cluster_F_specific_mouse[!is.na(cluster_F_specific_mouse)]

F_intersect <- intersect(opc12_uniqueF,cluster_F_specific_mouse)
F_opc12_only <- setdiff(opc12_uniqueF,cluster_F_specific_mouse)
F_cluster_only <- setdiff(cluster_F_specific_mouse,opc12_uniqueF)
F_none <- setdiff(unique(orthologs$MGI.symbol),union(opc12_uniqueF,cluster_F_specific_mouse))

fisher.test(x = matrix(data = c(length(F_intersect),length(F_opc12_only),
                                length(F_cluster_only),length(F_none)),
                       nrow = 2,ncol = 2))

phm_annotation <- data.frame(row.names = colnames(opc12_rheg_mat),
                             F_cluster = as.character(row.names(opc12_rheg_mat) %in% F_intersect))
