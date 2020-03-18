x <- -1*as.numeric(opc12_projection1[opc12$info$type == "ten-cell"])

hist(x,breaks = 24)


sapply(list.files("./external/Stochprof/R",full.names = T),source)
library(numDeriv)
library(MASS)

x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 30)
stochasticProfilingML()
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_PC mu_2_gene_PC        sigma 
# 0.4662       1.8490     -12.6330       0.2440 
# 
# Value of negative log-likelihood function at MLE:
#   214.4106 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   444.9226 
# 
# Approx. 95% confidence intervals for MLE:
#                   lower       upper
# p_1            0.3833369   0.5509712
# mu_1_gene_PC   1.6875435   2.0104565
# mu_2_gene_PC -12.6330020 -12.6329980
# sigma          0.1404649   0.4238498
# 
# Top parameter combinations:
#   p_1 mu_1_gene_PC mu_2_gene_PC sigma   target
# [1,] 0.4662        1.849      -12.633 0.244 214.4106
# [2,] 0.4661        1.849      -12.035 0.244 214.4106
# [3,] 0.4662        1.849      -12.044 0.244 214.4106
# [4,] 0.4661        1.849      -12.704 0.244 214.4106
# [5,] 0.4662        1.849      -11.914 0.244 214.4106
# [6,] 0.4662        1.849      -11.755 0.244 214.4106

x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 40)
stochasticProfilingML()
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.5572           1.9480          -1.6510           0.1680 
# 
# Value of negative log-likelihood function at MLE:
#   213.9693 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   444.04 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1               0.49488120 0.6177687
# mu_1_gene_Gene 1  1.80937901 2.0866210
# mu_2_gene_Gene 1 -4.85063200 1.5486320
# sigma             0.08665747 0.3256961
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.5572            1.948           -1.651 0.168 213.9693
# [2,] 0.5572            1.948           -1.652 0.168 213.9693
# [3,] 0.5572            1.948           -1.650 0.168 213.9693
# [4,] 0.5571            1.948           -1.652 0.168 213.9693
# [5,] 0.5574            1.948           -1.645 0.168 213.9693
# [6,] 0.5573            1.948           -1.644 0.168 213.9693

x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 50)
stochasticProfilingML()
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# 
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.6522           2.0110          -0.9530           0.1340 
# 
# Value of negative log-likelihood function at MLE:
#   213.3788 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   442.859 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1               0.60377291 0.6976720
# mu_1_gene_Gene 1  1.95342624 2.0685738
# mu_2_gene_Gene 1 -2.34070350 0.4347035
# sigma             0.08759039 0.2049997
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.6522            2.011           -0.953 0.134 213.3788
# [2,] 0.6522            2.011           -0.955 0.134 213.3788
# [3,] 0.6522            2.011           -0.954 0.134 213.3788
# [4,] 0.6521            2.011           -0.954 0.134 213.3788
# [5,] 0.6521            2.011           -0.953 0.134 213.3788
# [6,] 0.6522            2.011           -0.952 0.134 213.3788

x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 60)
stochasticProfilingML()
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.7571           2.0710         -16.3130           0.1230 
# 
# Value of negative log-likelihood function at MLE:
#   212.6898 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   441.481 
# 
# Approx. 95% confidence intervals for MLE:
#   lower       upper
# p_1                0.7133258   0.7961007
# mu_1_gene_Gene 1   2.0405389   2.1014611
# mu_2_gene_Gene 1 -16.3130000 -16.3130000
# sigma              0.0859805   0.1759585
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.7571            2.071          -16.313 0.123 212.6898
# [2,] 0.7571            2.071          -13.884 0.123 212.6898
# [3,] 0.7571            2.071          -16.388 0.123 212.6898
# [4,] 0.7571            2.071          -15.024 0.123 212.6898
# [5,] 0.7571            2.071          -15.550 0.123 212.6898
# [6,] 0.7571            2.071          -16.171 0.123 212.6898

x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 75)
stochasticProfilingML()
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# 
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.8446           2.1800          -1.7540           0.1200 
# 
# Value of negative log-likelihood function at MLE:
#   212.2329 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   440.5672 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1               0.80246642 0.8791005
# mu_1_gene_Gene 1  2.14932417 2.2106758
# mu_2_gene_Gene 1 -6.89765796 3.3896580
# sigma             0.08446417 0.1704865
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.8446             2.18           -1.754  0.12 212.2329
# [2,] 0.8445             2.18           -1.756  0.12 212.2329
# [3,] 0.8445             2.18           -1.757  0.12 212.2329
# [4,] 0.8446             2.18           -1.755  0.12 212.2329
# [5,] 0.8445             2.18           -1.752  0.12 212.2329
# [6,] 0.8445             2.18           -1.750  0.12 212.2329

x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 90)
stochasticProfilingML()
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.8449           2.3370           0.5050           0.1030 
# 
# Value of negative log-likelihood function at MLE:
#   212.2326 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   440.5666 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1               0.80310598 0.8791574
# mu_1_gene_Gene 1  2.31101430 2.3629857
# mu_2_gene_Gene 1 -0.04698774 1.0569877
# sigma             0.07137377 0.1486400
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.8449            2.337            0.505 0.103 212.2326
# [2,] 0.8448            2.337            0.506 0.103 212.2326
# [3,] 0.8449            2.337            0.506 0.103 212.2326
# [4,] 0.8449            2.337            0.504 0.103 212.2326
# [5,] 0.8448            2.337            0.505 0.103 212.2326
# [6,] 0.8449            2.337            0.505 0.103 212.2330

x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 120)
y <- stochasticProfilingML()
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.8452           2.5930           1.5350           0.0790 
# 
# Value of negative log-likelihood function at MLE:
#   212.2336 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   440.5686 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.80169854 0.8805800
# mu_1_gene_Gene 1 2.57343850 2.6125615
# mu_2_gene_Gene 1 1.34173518 1.7282648
# sigma            0.05483214 0.1138201
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.8452            2.593            1.535 0.079 212.2336
# [2,] 0.8453            2.593            1.535 0.079 212.2336
# [3,] 0.8453            2.593            1.536 0.079 212.2336
# [4,] 0.8452            2.593            1.536 0.079 212.2336
# [5,] 0.8453            2.592            1.535 0.079 212.2336
# [6,] 0.8452            2.592            1.536 0.079 212.2337



# Flip direction
x <- as.numeric(-1*opc12_projection1[opc12$info$type == "ten-cell"] + 120)
y <- stochasticProfilingML()
