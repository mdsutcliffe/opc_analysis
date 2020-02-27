
source("./opc12/RHEGs_opc12.R")
source("./opc90/RHEGs_opc90.R")


write.table(x = opc12_uniqueF,file = "./data/opc12_uniqueF.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc12_uniqueM,file = "./data/opc12_uniqueM.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc12_rheg,file = "./data/opc12_rheg.txt",quote = F,sep = "\n",row.names = F,col.names = F)

write.table(x = opc90_uniqueF,file = "./data/opc90_uniqueF.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc90_uniqueM,file = "./data/opc90_uniqueM.txt",quote = F,sep = "\n",row.names = F,col.names = F)
write.table(x = opc90_rheg,file = "./data/opc90_rheg.txt",quote = F,sep = "\n",row.names = F,col.names = F)

length(opc12_rheg) + length(opc12_uniqueF)
length(opc12_rheg) + length(opc12_uniqueM)

sqrt(length(opc12_rheg) + length(opc12_uniqueF)) / sqrt(length(opc12_rheg) + length(opc12_uniqueM))


areaOfIntersection <- function(x0,y0,r0,x1,y1,r1) {
  rr0 <- r0^2
  rr1 <- r1^2
  c =sqrt((x1-x0)^2 + (y1-y0)^2)
  phi <- acos((rr0+c^2-rr1) / (2*r0*c))*2
  theta <- acos((rr1+c^2-rr0) / (2*r1*c))*2
  area1 <- 0.5*theta*rr1 - 0.5*rr1*sin(theta)
  area2 <- 0.5*phi*rr0 - 0.5*rr0*sin(phi)
  return(area1 + area2)
}

areaOfIntersection(x0 = 1.222722/2,y0 = 1.222722/2,r0 = 1.222722/2,
                   x1 = 1.365,y1 = 0.5,r1 = 1/2) * (length(opc12_rheg) + length(opc12_uniqueM))
