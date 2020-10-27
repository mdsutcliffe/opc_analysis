source("./import/import_opcBulk.R")

# Trp53
cko12 <- as.numeric(bulk$log2[bulk$log2$symbol == "Trp53",9 + which(bulk$info$day == 12 & bulk$info$genotype == "CKO")])
wt12 <- as.numeric(bulk$log2[bulk$log2$symbol == "Trp53",9 + which(bulk$info$day == 12 & bulk$info$genotype == "WT")])

cko90 <- as.numeric(bulk$log2[bulk$log2$symbol == "Trp53",9 + which(bulk$info$day == 90 & bulk$info$genotype == "CKO")])
wt90 <- as.numeric(bulk$log2[bulk$log2$symbol == "Trp53",9 + which(bulk$info$day == 90 & bulk$info$genotype == "WT")])

cko150 <- as.numeric(bulk$log2[bulk$log2$symbol == "Trp53",9 + which(bulk$info$day == 150 & bulk$info$genotype == "CKO")])
wt150 <- as.numeric(bulk$log2[bulk$log2$symbol == "Trp53",9 + which(bulk$info$day == 150 & bulk$info$genotype == "WT")])

pdf(file = "./new/p53.pdf",width = 3.125,height = 2,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = NA,y = NA,
     xlim = c(0,9),ylim = c(0,10),
     las = 1,
     frame = F,axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "",ylab = "Log2(TPM + 1)",
     main = "Trp53")
points(x = rep(1,length(wt12)),
       y = wt12)
points(x = rep(2,length(cko12)),
       y = cko12)
points(x = rep(4,length(wt90)),
       y = wt90)
points(x = rep(5,length(cko90)),
       y = cko90)
points(x = rep(7,length(wt150)),
       y = wt150)
points(x = rep(8,length(cko150)),
       y = cko150)
axis(side = 2,lwd = 0.5/0.75,las = 1)
axis(side = 1,at = c(1,2,4,5,7,8),labels = c("WT","CKO","WT","CKO","WT","CKO"),las = 1,lwd = 0.5/0.75)
axis(side = 1,at = c(1.5,4.5,7.5),labels = c("12 dpi","90 dpi","150 dpi"),las = 1,lwd = 0.5/0.75,lwd.ticks = 0,padj = 2)
dev.off()


# Nf1
cko12 <- as.numeric(bulk$log2[bulk$log2$symbol == "Nf1",9 + which(bulk$info$day == 12 & bulk$info$genotype == "CKO")])
wt12 <- as.numeric(bulk$log2[bulk$log2$symbol == "Nf1",9 + which(bulk$info$day == 12 & bulk$info$genotype == "WT")])

cko90 <- as.numeric(bulk$log2[bulk$log2$symbol == "Nf1",9 + which(bulk$info$day == 90 & bulk$info$genotype == "CKO")])
wt90 <- as.numeric(bulk$log2[bulk$log2$symbol == "Nf1",9 + which(bulk$info$day == 90 & bulk$info$genotype == "WT")])

cko150 <- as.numeric(bulk$log2[bulk$log2$symbol == "Nf1",9 + which(bulk$info$day == 150 & bulk$info$genotype == "CKO")])
wt150 <- as.numeric(bulk$log2[bulk$log2$symbol == "Nf1",9 + which(bulk$info$day == 150 & bulk$info$genotype == "WT")])

pdf(file = "./new/Nf1.pdf",width = 3.125,height = 2,pointsize = 7,useDingbats = F)
par(mai = c(0.5,0.5,0.25,0.25),mgp = c(1.6,0.6,0),xpd = T,lwd = 0.5/0.75)
plot(x = NA,y = NA,
     xlim = c(0,9),ylim = c(0,10),
     las = 1,
     frame = F,axes = F,
     xaxs = "i",yaxs = "i",
     xlab = "",ylab = "Log2(TPM + 1)",
     main = "Nf1")
points(x = rep(1,length(wt12)),
       y = wt12)
points(x = rep(2,length(cko12)),
       y = cko12)
points(x = rep(4,length(wt90)),
       y = wt90)
points(x = rep(5,length(cko90)),
       y = cko90)
points(x = rep(7,length(wt150)),
       y = wt150)
points(x = rep(8,length(cko150)),
       y = cko150)
axis(side = 2,lwd = 0.5/0.75,las = 1)
axis(side = 1,at = c(1,2,4,5,7,8),labels = c("WT","CKO","WT","CKO","WT","CKO"),las = 1,lwd = 0.5/0.75)
axis(side = 1,at = c(1.5,4.5,7.5),labels = c("12 dpi","90 dpi","150 dpi"),las = 1,lwd = 0.5/0.75,lwd.ticks = 0,padj = 2)
dev.off()
