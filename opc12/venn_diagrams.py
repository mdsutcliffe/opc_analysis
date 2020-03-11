# Figure 3F

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

fig = plt.figure(figsize=(4,2))
plt.rcParams["font.family"] = "Arial"
v = venn2(subsets = (569,919,138),set_labels = ("Male","Female"),set_colors = ("#00000000","#00000000"))
plt.savefig("../plots/venn_opc12_sex.pdf",transparent = True)

fig = plt.figure(figsize=(4,2))
plt.rcParams["font.family"] = "Arial"
v = venn2(subsets = (134,465,4),set_labels = ("Stochastic profiling","Bulk RNA-seq"),set_colors = ("#00000000","#00000000"))
plt.savefig("../plots/venn_opc12_10c_vs_bulk.pdf",transparent = True)
