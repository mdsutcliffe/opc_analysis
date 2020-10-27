import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

os.chdir("/Users/mdsutcliffe/Github/opc_analysis")

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 10
plt.rcParams['pdf.fonttype'] = 42


fig = plt.figure(figsize=(3.33,2))
v = venn2(subsets = (137,651,63),set_labels = ("EMT Hallmark","PC Subtype"))
plt.savefig("./plots/venn_emt_pca.pdf",transparent = True)
