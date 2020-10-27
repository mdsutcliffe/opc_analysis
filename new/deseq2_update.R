# Sex-specific DE genes vs. combined model

load("./build/opcBulk_DESeq2_results.RData")

de12 <- bulk$deseq2$de12$genesDE
de90 <- bulk$deseq2$de90$genesDE

load("./new/opcBulk_DESeq2_results_revision.RData")

de12F <- bulk$deseq2$de12F$genesDE
de12M <- bulk$deseq2$de12M$genesDE
de90F <- bulk$deseq2$de90F$genesDE
de90M <- bulk$deseq2$de90M$genesDE
de150F <- bulk$deseq2$de150F$genesDE
de150M <- bulk$deseq2$de150M$genesDE

A <- de12
B <- de12F
C <- de12M
length(setdiff(A,union(B,C)))
length(setdiff(B,union(A,C)))
length(setdiff(intersect(A,B),C))
length(setdiff(C,union(A,B)))
length(setdiff(intersect(A,C),B))
length(setdiff(intersect(B,C),A))
length(intersect(A,intersect(B,C)))

cat(intersect(A,intersect(B,C)))

A <- de90
B <- de90F
C <- de90M
length(setdiff(A,union(B,C)))
length(setdiff(B,union(A,C)))
length(setdiff(intersect(A,B),C))
length(setdiff(C,union(A,B)))
length(setdiff(intersect(A,C),B))
length(setdiff(intersect(B,C),A))
length(intersect(A,intersect(B,C)))

cat(intersect(A,intersect(B,C)))

# Sex-specific DE genes overlap with RHEGs

source("./functions/opc12_rheg.R")

A <- union(opc12$candidates$Fspecific,opc12$candidates$RHEG)
B <- de12F
length(setdiff(A,B))
length(setdiff(B,A))
length(intersect(A,B))

cat(intersect(A,B))

A <- union(opc12$candidates$Mspecific,opc12$candidates$RHEG)
B <- de12M
length(setdiff(A,B))
length(setdiff(B,A))
length(intersect(A,B))

cat(intersect(A,B))

source("./functions/opc90_rheg.R")

A <- union(opc90$candidates$Fspecific,opc90$candidates$RHEG)
B <- de90F
length(setdiff(A,B))
length(setdiff(B,A))
length(intersect(A,B))

cat(intersect(A,B))

A <- union(opc90$candidates$Mspecific,opc90$candidates$RHEG)
B <- de90M
length(setdiff(A,B))
length(setdiff(B,A))
length(intersect(A,B))

cat(intersect(A,B))