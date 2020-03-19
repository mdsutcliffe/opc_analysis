# Figure 4D
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

os.chdir("/Users/mdsutcliffe/Github/opc_analysis")

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 10
plt.rcParams['pdf.fonttype'] = 42

# Left, male-female overlap
fig = plt.figure(figsize=(2,2))
v = venn2(subsets = (619,1203,231),set_labels = ("Male","Female"))
plt.savefig("./plots/venn_opc90_sex.pdf",transparent = True)

# Right
fig = plt.figure(figsize=(2,2))
v = venn2(subsets = (224,136,7),set_labels = ("Stochastic profiling","Bulk RNA-seq"))
plt.savefig("./plots/venn_opc90_10c_vs_bulk.pdf",transparent = True)
