# Figure 3F
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

os.chdir("/Users/mdsutcliffe/Github/opc_analysis")

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 10
plt.rcParams['pdf.fonttype'] = 42

# Left, male-female overlap
fig = plt.figure(figsize=(2,2))
v = venn2(subsets = (569,919,138),set_labels = ("Male","Female"))
plt.savefig("./plots/venn_opc12_sex.pdf",transparent = True)

# Right, stochastic profiling-bulk rna-seq overlap
fig = plt.figure(figsize=(2,2))
v = venn2(subsets = (134,465,4),set_labels = ("Stochastic profiling","Bulk RNA-seq"))
plt.savefig("./plots/venn_opc12_10c_vs_bulk.pdf",transparent = True)
