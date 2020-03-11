# Figure 4D

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

fig = plt.figure(figsize=(4,2))
plt.rcParams["font.family"] = "Arial"
v = venn2(subsets = (619,1203,231),set_labels = ("Male","Female"),set_colors = ("#00000000","#00000000"))
plt.savefig("../plots/venn_opc90_sex.pdf",transparent = True)

fig = plt.figure(figsize=(4,2))
plt.rcParams["font.family"] = "Arial"
v = venn2(subsets = (224,136,7),set_labels = ("Stochastic profiling","Bulk RNA-seq"),set_colors = ("#00000000","#00000000"))
plt.savefig("../plots/venn_opc90_10c_vs_bulk.pdf",transparent = True)
