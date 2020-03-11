x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 30)

hist(x,breaks = 24)


sapply(list.files("./external/Stochprof/R",full.names = T),source)
library(numDeriv)
library(MASS)

x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 30)
stochasticProfilingML()
# ***** Estimation started! *****
#   
#   
#   Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.1556           2.1500           0.7650           0.4110 
# 
# Value of negative log-likelihood function at MLE:
#   203.2877 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   422.6768 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.04583754 0.4141235
# mu_1_gene_Gene 1 1.59899310 2.7010069
# mu_2_gene_Gene 1 0.51858244 1.0114176
# sigma            0.23002435 0.7343614
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.1556            2.150            0.765 0.411 203.2877
# [2,] 0.1553            2.151            0.765 0.411 203.2877
# [3,] 0.1554            2.151            0.765 0.411 203.2877
# [4,] 0.1555            2.150            0.765 0.411 203.2877
# [5,] 0.1552            2.151            0.765 0.411 203.2877
# [6,] 0.1554            2.150            0.765 0.411 203.2877

x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 40)
stochasticProfilingML()
# ***** Estimation started! *****
#   
#   
#   Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.1198           2.4160           1.1880           0.3020 
# 
# Value of negative log-likelihood function at MLE:
#   203.092 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   422.2854 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.05206391 0.2522146
# mu_1_gene_Gene 1 2.02517603 2.8068240
# mu_2_gene_Gene 1 1.08376006 1.2922399
# sigma            0.17722718 0.5146163
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma  target
# [1,] 0.1198            2.416            1.188 0.302 203.092
# [2,] 0.1197            2.416            1.188 0.302 203.092
# [3,] 0.1197            2.416            1.189 0.302 203.092
# [4,] 0.1199            2.415            1.188 0.302 203.092
# [5,] 0.1199            2.416            1.188 0.301 203.092
# [6,] 0.1198            2.415            1.188 0.302 203.092

x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 50)
stochasticProfilingML()
# ***** Estimation started! *****
#   
#   
#   Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.1091           2.5650           1.4680           0.2300 
# 
# Value of negative log-likelihood function at MLE:
#   202.9478 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   421.997 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.06060102 0.1886194
# mu_1_gene_Gene 1 2.30132224 2.8286778
# mu_2_gene_Gene 1 1.40947100 1.5265290
# sigma            0.13707462 0.3859212
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.1091            2.565            1.468 0.230 202.9478
# [2,] 0.1090            2.565            1.468 0.230 202.9478
# [3,] 0.1092            2.565            1.468 0.230 202.9478
# [4,] 0.1092            2.564            1.468 0.230 202.9478
# [5,] 0.1112            2.576            1.463 0.222 202.9845
# [6,] 0.1019            2.599            1.471 0.257 203.1108

x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 60)
stochasticProfilingML()
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.1052           2.6660           1.6810           0.1840 
# 
# Value of negative log-likelihood function at MLE:
#   202.8201 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   421.7416 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.06538205 0.1649866
# mu_1_gene_Gene 1 2.47270487 2.8592951
# mu_2_gene_Gene 1 1.64056130 1.7214387
# sigma            0.11423223 0.2963787
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.1052            2.666            1.681 0.184 202.8201
# [2,] 0.1051            2.666            1.681 0.184 202.8201
# [3,] 0.1051            2.667            1.681 0.184 202.8201
# [4,] 0.1053            2.666            1.681 0.184 202.8201
# [5,] 0.0982            2.695            1.684 0.190 202.8811
# [6,] 0.1040            2.659            1.684 0.174 202.8942

x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 75)
stochasticProfilingML()
# ***** Estimation started! *****
#   
#   
#   Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.1029           2.7810           1.9310           0.1420 
# 
# Value of negative log-likelihood function at MLE:
#   202.6681 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   421.4376 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.06817474 0.1524196
# mu_1_gene_Gene 1 2.63771053 2.9242895
# mu_2_gene_Gene 1 1.90250301 1.9594970
# sigma            0.09290989 0.2170275
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.1029            2.781            1.931 0.142 202.6681
# [2,] 0.1030            2.781            1.931 0.142 202.6681
# [3,] 0.1031            2.781            1.931 0.142 202.6682
# [4,] 0.1028            2.782            1.931 0.142 202.6682
# [5,] 0.1029            2.781            1.931 0.142 202.6683
# [6,] 0.1034            2.780            1.932 0.145 202.6809

x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 90)
stochasticProfilingML()
# ***** Estimation started! *****
#   
#   [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# 
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.1021           2.8760           2.1290           0.1170 
# 
# Value of negative log-likelihood function at MLE:
#   202.5648 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   421.231 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.06897823 0.1485878
# mu_1_gene_Gene 1 2.75491572 2.9970843
# mu_2_gene_Gene 1 2.10612015 2.1518798
# sigma            0.07826159 0.1749134
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.1021            2.876            2.129 0.117 202.5648
# [2,] 0.1020            2.876            2.129 0.117 202.5648
# [3,] 0.1022            2.876            2.129 0.117 202.5648
# [4,] 0.1021            2.876            2.129 0.117 202.5653
# [5,] 0.0984            2.889            2.128 0.118 202.6041
# [6,] 0.0988            2.878            2.133 0.125 202.6379

x <- as.numeric(opc90_projection1[opc90$info$type == "ten-cell"] + 120)
y <- stochasticProfilingML()
# ***** Estimation started! *****
#   
#   [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# [1] "Hessian not positive semi-definite; return NULL."
# 
# Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1            sigma 
# 0.1015           3.0360           2.4360           0.0870 
# 
# Value of negative log-likelihood function at MLE:
#   202.4496 
# 
# Violation of constraints:
#   none
# 
# BIC:
#   421.0006 
# 
# Approx. 95% confidence intervals for MLE:
#   lower     upper
# p_1              0.06890331 0.1470815
# mu_1_gene_Gene 1 2.93832726 3.1336727
# mu_2_gene_Gene 1 2.41899458 2.4530054
# sigma            0.05892816 0.1284445
# 
# Top parameter combinations:
#   p_1 mu_1_gene_Gene 1 mu_2_gene_Gene 1 sigma   target
# [1,] 0.1015            3.036            2.436 0.087 202.4496
# [2,] 0.1016            3.036            2.436 0.087 202.4496
# [3,] 0.1015            3.036            2.436 0.087 202.4497
# [4,] 0.1014            3.036            2.436 0.087 202.4497
# [5,] 0.1016            3.036            2.436 0.086 202.4497
# [6,] 0.1015            3.036            2.436 0.087 202.4500