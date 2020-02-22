# Check PC genes for chromosome

f <- "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Sham/janes_dataset2_rsem_counts.csv"

x <- read.csv(f,header = T)
x2 <- x

x3 <- x2[x2$symbol %in% pc$Gene,1:9]
x4 <- table(x3$chr)
x5 <- x4[nchar(names(x4)) <= 3]

y <- table(x$chr)
y2 <- y[nchar(names(y)) <= 3]

z <- x5/y2
z <- z[c(1,12,16:22,2:11,13:15,23:25)]
barplot(z,xlab = "chromosome",ylab = "# PC genes / Chromosome size")
