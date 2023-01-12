import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

base_path = Path(__file__).parent
outdir = '../img/'
kashif_version = {'nonorm': '../csv/kashif/nonorm/mode0/', 'norm': '../csv/kashif/norm/mode0/'}
curr_version = kashif_version['nonorm']
csv_files = {"kashif": (base_path / f"{curr_version}kashif_querytime.csv").resolve(), "pexeso": (base_path / "../csv/pexeso/tau6%-T1%/pexeso_querytime.csv").resolve()}

# plt.rcParams["figure.figsize"] = (8,7)
k = pd.read_csv(csv_files['kashif'])
p = pd.read_csv(csv_files['pexeso'])


dfk = k.groupby(
    ['k', 'nb_threads']
).agg(
    querytime = ('querytime','mean'),
).reset_index()

dfp = p.groupby(
    ['tau', 'join_threashold']
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
plt.subplots_adjust(top=0.9)
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(left=0.15)
plt.title("Mean querytime Kashif vs PEXESO (100k tables, 4.9M vectors)")
plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
plt.ylabel(r'$\mathrm{mean\ querytime\ (sec)}$', fontsize = 11)
# plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.yscale("log")
plt.grid()  #just add this
plt.savefig(f"{outdir}mean_querytime_kashif_vs_pexeso(nonorm,log).png")
plt.close()
