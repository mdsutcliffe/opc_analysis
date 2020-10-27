setwd('/Volumes/JanesLab/Matthew Sutcliffe/Image analysis')
list.files()

froot1 <- 'F6340_M8170-2_gH2AX-Upf3b_'
x1 <- read.csv(file = paste0(froot1,'opc_expand.csv'))
xinfo1 <- read.csv(file = paste0(froot1,'Image.csv'))

froot2 <- 'F7460_M8170-1_Upf3b-gH2AX_'
x2 <- read.csv(file = paste0(froot2,'opc_expand.csv'))
xinfo2 <- read.csv(file = paste0(froot2,'Image.csv'))

hist(x1$Intensity_MedianIntensity_c488_rescale_smooth[grepl('F6340',xinfo1$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "gH2AX_median_intensity",ylab = "Number of cells",main = "F6340",las = 1)
hist(x1$Intensity_MedianIntensity_c488_rescale_smooth[grepl('M8170',xinfo1$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "gH2AX_median_intensity",ylab = "Number of cells",main = "M8170_OB#2",las = 1)
hist(x2$Intensity_MedianIntensity_c647_rescale_smooth[grepl('F7460',xinfo2$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "gH2AX_median_intensity",ylab = "Number of cells",main = "F7460",las = 1)
hist(x2$Intensity_MedianIntensity_c647_rescale_smooth[grepl('M8170',xinfo2$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "gH2AX_median_intensity",ylab = "Number of cells",main = "M8170_OB#1",las = 1)

hist(x1$Intensity_MedianIntensity_c647_rescale_smooth[grepl('F6340',xinfo1$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "Upf3b_median_intensity",ylab = "Number of cells",main = "F6340",las = 1)
hist(x1$Intensity_MedianIntensity_c647_rescale_smooth[grepl('M8170',xinfo1$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "Upf3b_median_intensity",ylab = "Number of cells",main = "M8170_OB#2",las = 1)
hist(x2$Intensity_MedianIntensity_c488_rescale_smooth[grepl('F7460',xinfo2$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "Upf3b_median_intensity",ylab = "Number of cells",main = "F7460",las = 1)
hist(x2$Intensity_MedianIntensity_c488_rescale_smooth[grepl('M8170',xinfo2$FileName_c488)],breaks = seq(0,1,0.01),
     xlab = "Upf3b_median_intensity",ylab = "Number of cells",main = "M8170_OB#1",las = 1)

plot(x = x1$Intensity_MedianIntensity_c488_rescale_smooth[grepl('F6340',xinfo1$FileName_c488)],
     y = x1$Intensity_MedianIntensity_c647_rescale_smooth[grepl('F6340',xinfo1$FileName_c488)],
     xlim = c(0,0.5),ylim = c(0,0.5),xlab = "gH2AX_median_intensity",ylab = "Upf3b_median_intensity",main = "F6340",las = 1)
plot(x = x1$Intensity_MedianIntensity_c488_rescale_smooth[grepl('M8170',xinfo1$FileName_c488)],
     y = x1$Intensity_MedianIntensity_c647_rescale_smooth[grepl('M8170',xinfo1$FileName_c488)],
     xlim = c(0,0.5),ylim = c(0,0.5),xlab = "gH2AX_median_intensity",ylab = "Upf3b_median_intensity",main = "M8170_OB#2",las = 1)


plot(x = x2$Intensity_MedianIntensity_c647_rescale_smooth[grepl('F7460',xinfo2$FileName_c488)],
     y = x2$Intensity_MedianIntensity_c488_rescale_smooth[grepl('F7460',xinfo2$FileName_c488)],
     xlim = c(0,0.5),ylim = c(0,0.5),xlab = "gH2AX_median_intensity",ylab = "Upf3b_median_intensity",main = "F7460",las = 1)
plot(x = x2$Intensity_MedianIntensity_c647_rescale_smooth[grepl('M8170',xinfo2$FileName_c488)],
     y = x2$Intensity_MedianIntensity_c488_rescale_smooth[grepl('M8170',xinfo2$FileName_c488)],
     xlim = c(0,0.5),ylim = c(0,0.5),xlab = "gH2AX_median_intensity",ylab = "Upf3b_median_intensity",main = "M8170_OB#1",las = 1)
