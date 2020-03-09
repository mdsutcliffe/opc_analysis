x <- as.numeric(opc12_projection1[opc12$info$type == "ten-cell"] + 60)

hist(x,breaks = seq(30,90,3.3333))


# ***** Estimation started! *****
#   
#   
#   Maximum likelihood estimate (MLE):
#   p_1 mu_1_gene_PC1 mu_2_gene_PC1         sigma 
# 0.7571        2.0710       -9.9380        0.1230 
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
#   lower      upper
# p_1            0.71332532  0.7961010
# mu_1_gene_PC1  2.04053845  2.1014616
# mu_2_gene_PC1 -9.93800138 -9.9379986
# sigma          0.08598063  0.1759582
# 
# Top parameter combinations:
#   p_1 mu_1_gene_PC1 mu_2_gene_PC1 sigma   target
# [1,] 0.7571         2.071        -9.938 0.123 212.6898
# [2,] 0.7571         2.071        -9.266 0.123 212.6898
# [3,] 0.7571         2.071        -7.140 0.123 212.6898
# [4,] 0.7571         2.071        -9.097 0.123 212.6898
# [5,] 0.7571         2.071        -9.124 0.123 212.6898
# [6,] 0.7571         2.071       -10.465 0.123 212.6898



x <- exp(as.numeric(opc12_projection1[opc12$info$type == "ten-cell"]) / 20)

stochasticProfilingML()
