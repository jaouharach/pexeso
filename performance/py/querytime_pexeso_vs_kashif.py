import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

base_path = Path(__file__).parent
# outdir = '../img/'
outdir = '../img_(k=10k)/'
k = 'k=10k'

kashif_versions = {'nonorm': f'../csv/kashif/{k}/nonorm/mode0/', 'norm': f'../csv/kashif/{k}/norm/mode0/'}
curr_version = kashif_versions['nonorm']
pexeso_versions =  {(0, 'violet'): 'tau0%-T1%/', (2, 'magenta'): 'tau2%-T1%/', (6, 'indigo'): 'tau6%-T1%/'}
csv_files = {"kashif": (base_path / f"{curr_version}kashif_querytime.csv").resolve(), "pexeso": (base_path / f"../csv/pexeso/").resolve()}

# plot kashif query time
k = pd.read_csv(csv_files['kashif'])
dfk = k.groupby(
    ['k', 'nb_threads']
).agg(
    querytime = ('querytime','mean'),
).reset_index()
print(dfk)
g = sns.catplot(data=dfk, x='k', y='querytime', hue='nb_threads',  kind='point', scale = 0.5,
    palette=sns.color_palette('nipy_spectral_r', 3),  markers=['o', 'v', '*', 'X'], legend=False)

# plot pexeso query time
for (v, c) in pexeso_versions:
    csv_file = str(csv_files['pexeso']) + f"/{pexeso_versions[(v,c)]}pexeso_querytime.csv"
    p = pd.read_csv(csv_file)
    dfp = p.groupby(
        ['tau', 'join_threashold']
    ).agg(
        querytime = ('querytime','mean'),
    ).reset_index()
    print(dfp)

    g.map(plt.axhline, y=dfp['querytime'][0], ls='--', c=f'{c}', label=f'pexeso (tau = {v}%)')
    g.set_xticklabels(rotation=30)

plt.subplots_adjust(top=0.9)
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(left=0.15)
plt.title("Mean querytime Kashif (k = 10.000) vs PEXESO \n(100k tables, 4.9M vectors)")
plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
plt.ylabel(r'$\mathrm{mean\ querytime\ (sec)}$', fontsize = 11)
plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.yscale("log")
plt.grid()  #just add this
plt.savefig(f"{outdir}mean_querytime_kashif_vs_pexeso(nonorm,log).png")
plt.close()
