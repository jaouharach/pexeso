import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

base_path = Path(__file__).parent
outdir = '../img_(k=1k)/'
k = 'k=1k'
kashif_versions = {'nonorm': f'../csv/kashif/{k}/nonorm/', 'norm': f'../csv/kashif/{k}/norm/'}
search_modes = ["mode0/"]
norm_csv = (base_path / f"{kashif_versions['norm']}{search_modes[0]}/kashif_querytime.csv").resolve()
nonorm_csv = (base_path / f"{kashif_versions['nonorm']}{search_modes[0]}/kashif_querytime.csv").resolve()

k_norm = pd.read_csv(norm_csv)
k_nonorm = pd.read_csv(nonorm_csv)


dfk = k_nonorm.groupby(
    ['k', 'nb_threads']
).agg(
    nonorm = ('querytime','mean'),
).reset_index()

dfk_norm = k_norm.groupby(
    ['k', 'nb_threads']
).agg(
    norm = ('querytime','mean'),
).reset_index()

dfk["norm"] = dfk_norm["norm"] 

dfk.drop(dfk[dfk['nb_threads'] != 64].index, inplace = True)
dfm = dfk.melt(['k', 'nb_threads'], var_name='data', value_name='mean_querytime')

print(dfm)

plt.rcParams["figure.figsize"] = (15,13)
flatui = ["#e74c3c", "#34495e"]
g = sns.catplot(data=dfm, x='k', y='mean_querytime', hue='data',  kind='point', scale = 0.5,
    palette=sns.color_palette(flatui),  markers=['o', 'v'])

plt.subplots_adjust(top=0.9)
# plt.subplots_adjust(left=0.15)
plt.title("Kashif mean querytime (100k tables, 4.9M vectors, 64 threads)")
plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
plt.ylabel(r'$\mathrm{mean\ querytime\ (sec)}$', fontsize = 11)
# plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"{outdir}kashif_mean_querytime(norm vs no norm).png")
plt.close()
