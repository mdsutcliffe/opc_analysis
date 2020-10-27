import os
import matplotlib.pyplot as plt
import matplotlib_venn as venn

os.chdir("/Users/mdsutcliffe/Github/opc_analysis")

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 10
plt.rcParams['pdf.fonttype'] = 42


fig = plt.figure(figsize=(4,3))
v = venn.venn3(subsets = (278,98,110,57,66,0,15),set_labels = ("DE12_factor","DE12_Female","DE12_Male"))
plt.savefig("./plots/DE_sex_vs_factor_12.pdf",transparent = True)

fig = plt.figure(figsize=(4,3))
v = venn.venn3(subsets = (78,15,14,156,40,2,11),set_labels = ("DE90_factor","DE90_Female","DE90_Male"))
plt.savefig("./plots/DE_sex_vs_factor_90.pdf",transparent = True)

fig = plt.figure(figsize=(6,3))
v = venn.venn2(subsets = (1041,207,16),set_labels = ("Female_heterogeneities_12","DE12_Female"))
plt.savefig("./plots/DE_heterogeneous_12F.pdf",transparent = True)
fig = plt.figure(figsize=(6,3))
v = venn.venn2(subsets = (700,131,7),set_labels = ("Male_heterogeneities_12","DE12_Male"))
plt.savefig("./plots/DE_heterogeneous_12M.pdf",transparent = True)

fig = plt.figure(figsize=(6,3))
v = venn.venn2(subsets = (1428,36,6),set_labels = ("Female_heterogeneities_90","DE90_Female"))
plt.savefig("./plots/DE_heterogeneous_90F.pdf",transparent = True)
fig = plt.figure(figsize=(6,3))
v = venn.venn2(subsets = (838,197,12),set_labels = ("Male_heterogeneities_90","DE90_Male"))
plt.savefig("./plots/DE_heterogeneous_90M.pdf",transparent = True)
