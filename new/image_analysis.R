f = "/Volumes/GoogleDrive/My Drive/Janes Lab/opc_filter.csv"
x = read.csv(file = f)

hist(x$Intensity_UpperQuartileIntensity_Actn1_rescale_smooth_backgroundSubtract,breaks = seq(0,1,0.01),xlim = c(0,0.5))
hist(x$Intensity_UpperQuartileIntensity_Fn1_rescale_smooth_backgroundSubtract,breaks = seq(0,1,0.01),xlim = c(0,0.5))

