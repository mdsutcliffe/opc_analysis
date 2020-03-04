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

# Human liberal
{
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
}

# Mouse liberal
{
cluster_F_specific_mouse <- orthologs$MGI.symbol[match(x = cluster_F_specific,table = orthologs$HGNC.symbol)]
cluster_F_specific_mouse <- cluster_F_specific_mouse[!is.na(cluster_F_specific_mouse)]

F_intersect <- intersect(opc12_uniqueF,cluster_F_specific_mouse)
F_opc12_only <- setdiff(opc12_uniqueF,cluster_F_specific_mouse)
F_cluster_only <- setdiff(cluster_F_specific_mouse,opc12_uniqueF)
F_none <- setdiff(unique(orthologs$MGI.symbol),union(opc12_uniqueF,cluster_F_specific_mouse))

fisher.test(x = matrix(data = c(length(F_intersect),length(F_opc12_only),
                                length(F_cluster_only),length(F_none)),
                       nrow = 2,ncol = 2))
}



phm_annotation <- data.frame(row.names = colnames(opc12_rheg_mat),
                             F_cluster = as.character(row.names(opc12_rheg_mat) %in% F_intersect))


# Mouse conservative
{
  cluster_F_specific_mouse <- orthologs$MGI.symbol[match(x = cluster_F_specific,table = orthologs$HGNC.symbol)]
  cluster_F_specific_mouse <- cluster_F_specific_mouse[!is.na(cluster_F_specific_mouse)]
  
  F_intersect <- intersect(opc12_uniqueF,cluster_F_specific_mouse)
  F_opc12_only <- setdiff(opc12_uniqueF,cluster_F_specific_mouse)
  F_cluster_only <- setdiff(cluster_F_specific_mouse,opc12_uniqueF)
  F_none <- setdiff(orthologs$MGI.symbol[!duplicated(orthologs$MGI.symbol) & !duplicated(orthologs$HGNC.symbol)],
                    union(opc12_uniqueF,cluster_F_specific_mouse))
  
  fisher.test(x = matrix(data = c(length(F_intersect),length(F_opc12_only),
                                  length(F_cluster_only),length(F_none)),
                         nrow = 2,ncol = 2))
  
  
  
  cluster_M_specific_mouse <- orthologs$MGI.symbol[match(x = cluster_M_specific,table = orthologs$HGNC.symbol)]
  cluster_M_specific_mouse <- cluster_M_specific_mouse[!is.na(cluster_M_specific_mouse)]
  
  M_intersect <- intersect(opc12_uniqueM,cluster_M_specific_mouse)
  M_opc12_only <- setdiff(opc12_uniqueM,cluster_M_specific_mouse)
  M_cluster_only <- setdiff(cluster_M_specific_mouse,opc12_uniqueM)
  M_none <- setdiff(orthologs$MGI.symbol[!duplicated(orthologs$MGI.symbol) & !duplicated(orthologs$HGNC.symbol)],
                    union(opc12_uniqueM,cluster_M_specific_mouse))
  
  fisher.test(x = matrix(data = c(length(M_intersect),length(M_opc12_only),
                                  length(M_cluster_only),length(M_none)),
                         nrow = 2,ncol = 2))
}


opc12_F_mat <- opc12$log2[match(x = opc12_uniqueF,table = opc12$log2$symbol),c(1:9,9+which(opc12$info$type == "ten-cell" & opc12$info$sex == "female"))]
row.names(opc12_F_mat) <- opc12_F_mat$symbol
opc12_F_mat <- opc12_F_mat[,10:ncol(opc12_F_mat)]
opc12_F_mat <- as.matrix(opc12_F_mat)

annotation_samples <- data.frame(row.names = colnames(opc12_F_mat),
                                   sex = opc12$info$sex[opc12$info$type == "ten-cell" & opc12$info$sex == "female"])

annotation_gene <- data.frame(row.names = row.names(opc12_F_mat),
                              F.cluster = as.character(row.names(opc12_F_mat) %in% cluster_F_specific_mouse))

annotation_colors <- list(F.cluster = c("TRUE" = "#000000","FALSE" = "#00000000"))

pdf(file = "./plots/RHEG_heatmap_opc12_sex_cluster_F.pdf",width = 8,height = 8)
par(mar=c(3.5,3.5,1,1),mgp = c(2.5,1,0))
pheatmap(mat = opc12_F_mat,
         color = rev(brewer.pal(11,"RdBu")),
         breaks = seq(-4,4,length.out = 12),
         border_color = NA,
         clustering_method = "ward.D2",
         show_rownames = F,
         annotation_row = annotation_gene,
         scale = "row",
         labels_col = rep(x = "    ",ncol(opc12_F_mat)),
         annotation_colors = annotation_colors)
grid.text(label = "Female 10-cell samples",y = 0.04)
grid.text(label = paste(nrow(opc12_F_mat),"F-specific heterogeneities"),x = 0.86,y = 0.4,rot = 90)
dev.off()




opc12_M_mat <- opc12$log2[match(x = opc12_uniqueM,table = opc12$log2$symbol),c(1:9,9+which(opc12$info$type == "ten-cell" & opc12$info$sex == "male"))]
row.names(opc12_M_mat) <- opc12_M_mat$symbol
opc12_M_mat <- opc12_M_mat[,10:ncol(opc12_M_mat)]
opc12_M_mat <- as.matrix(opc12_M_mat)



annotation_gene <- data.frame(row.names = row.names(opc12_M_mat),
                              M.cluster = as.character(row.names(opc12_M_mat) %in% cluster_M_specific_mouse))

annotation_colors <- list(M.cluster = c("TRUE" = "#000000","FALSE" = "#00000000"))

pdf(file = "./plots/RHEG_heatmap_opc12_sex_cluster_M.pdf",width = 8,height = 8)
par(mar=c(3.5,3.5,1,1),mgp = c(2.5,1,0))
pheatmap(mat = opc12_M_mat,
         color = rev(brewer.pal(11,"RdBu")),
         breaks = seq(-4,4,length.out = 12),
         border_color = NA,
         clustering_method = "ward.D2",
         show_rownames = F,
         annotation_row = annotation_gene,
         scale = "row",
         labels_col = rep(x = "    ",ncol(opc12_M_mat)),
         annotation_colors = annotation_colors)
grid.text(label = "Male 10-cell samples",y = 0.04)
grid.text(label = paste(nrow(opc12_M_mat),"M-specific heterogeneities"),x = 0.86,y = 0.4,rot = 90)
dev.off()