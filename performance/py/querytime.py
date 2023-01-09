import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

base_path = Path(__file__).parent
outdir = '../img/'
kashif_version = {'nonorm': '../csv/kashif/nonorm/', 'norm': '../csv/kashif/norm/'}
curr_version = kashif_version['norm']
csv_files = {"kashif": (base_path / f"{curr_version}kashif_querytime.csv").resolve(), "pexeso": (base_path / "../csv/pexeso/pexeso_querytime.csv").resolve()}

plt.rcParams["figure.figsize"] = (8,7)
k = pd.read_csv(csv_files['kashif'])
p = pd.read_csv(csv_files['pexeso'])


dfk = k.groupby(
    ['k', 'nb_threads']
).agg(
    querytime = ('querytime','mean'),
).reset_index()

dfp = p.groupby(
    ['tau', 'join_threashold', 'total_tables']
).agg(
    querytime = ('querytime','mean'),
).reset_index()

print(dfk)
print(dfp)

g = sns.catplot(data=dfk, x='k', y='querytime', hue='nb_threads',  kind='point', scale = 0.5,
    palette=sns.color_palette('nipy_spectral_r', 4),  markers=['o', 'v', '*', 'X'])

for t in dfp['querytime']:
    g.map(plt.axhline, y=t, ls='--', c='r')
    break

g.set_xticklabels(rotation=30)
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(left=0.15)
plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
plt.ylabel(r'$\mathrm{mean\ querytime\ (sec)}$', fontsize = 11)
# plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"{outdir}mean_querytime(norm).png")
plt.close()
