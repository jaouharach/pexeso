import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

base_path = Path(__file__).parent
outdir = '../img/'

csv_file = (base_path / f"../csv/pexeso/pexeso_querytime_all.csv").resolve()
p = pd.read_csv(csv_file)

dfp = p.groupby(
    ['tau', 'join_threashold']
).agg(
    mean_querytime = ('querytime','mean'),
).reset_index()

dfp['tau'] = dfp['tau'] * 100
dfp['tau'] = dfp['tau'].astype(str)

flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
# plt.rcParams["figure.figsize"] = (15,13)
g = sns.lineplot(data=dfp, x='tau', y='mean_querytime', color="#2ecc71",  marker='v')

# plt.subplots_adjust(top=0.9)
# plt.subplots_adjust(left=0.15)
plt.title("Pexeso mean querytime (100k tables, 4.9M vectors)")
plt.xlabel(r'$\mathrm{tau\ (\%)}$', fontsize = 11)
plt.ylabel(r'$\mathrm{mean\ querytime\ (sec)}$', fontsize = 11)
# plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"{outdir}pexeso_mean_querytime.png")
plt.close()
