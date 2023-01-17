import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import csv

# search_modes:
# 0 : "kashif returns all nn results"
# 1 : "kashif only returns NNs with the same distance to the query (as the first nn)"
# 2 : "similar to mode 1 but kashif also returns extra results found in the last increment"

base_path = Path(__file__).parent
outdir = '../img_(k=10k)/'
k = 'k=10k'
kashif_versions = {'nonorm': f'../csv/kashif/{k}/nonorm/', 'norm': f'../csv/kashif/{k}/norm/'}
pexeso_versions = {'0': '../csv/pexeso/tau0%-T1%/', '2': '../csv/pexeso/tau2%-T1%/', '6': '../csv/pexeso/tau6%-T1%/', '8': '../csv/pexeso/tau8%-T1%/'}

curr_kashif_version = kashif_versions['nonorm']
search_modes = ["mode0/", "mode1/", "mode2/"]

rows = []
for pv in pexeso_versions:
    curr_pexeso_version = pexeso_versions[pv]
    for sm in search_modes:
        csv_files = {"kashif": (base_path / f"{curr_kashif_version}{sm}kashif_results.csv").resolve(), "pexeso": (base_path / f"{curr_pexeso_version}pexeso_results.csv").resolve()}

        # plt.rcParams["figure.figsize"] = (8,7)
        k = pd.read_csv(csv_files['kashif'])
        p = pd.read_csv(csv_files['pexeso'])
        k = k.drop('overlap', axis=1)
        p = p.drop('nooverlap', axis=1)

        queryids = pd.unique(k['tqq'])
        # print(queryids)

        mean_recall = 0.0
        mean_precision = 0.0
        nb_query_results = 0

        for query in queryids:
            print(f"query: {query}")
            res = k.loc[k['tqq'] == query].reset_index(drop=True) # kashif results
            gt = p.loc[p['tqq'] == query].reset_index(drop=True) # ground truth results (pexeso)
            res['tp'] = res['tss'].map(gt['tss'].value_counts())
            res['tp'] = res['tp'].fillna(0)

            if gt.shape[0] != 0:
                recall = res['tp'].sum() / (gt.shape[0])
                precision = res['tp'].sum() / (res.shape[0])
                print(f"recall = {recall}, precision = {precision}")
                mean_recall += recall
                mean_precision += precision
                nb_query_results += 1

        mean_recall /= nb_query_results
        mean_precision /= nb_query_results
        rows.append({'tau' : pv, 'search-mode' : sm, 'avg. precision': mean_precision, 'avg. recall': mean_recall})


df = pd.DataFrame(rows)
print(df)

flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
g = sns.barplot(data=df, x='tau', y='avg. recall', hue='search-mode', palette=sns.color_palette(flatui))

plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(left=0.15)
plt.title("Kashif avg. recall (100k tables, 4.9M vectors, k = 10.000)")
plt.xlabel(r'$\mathrm{tau\ (\%)}$', fontsize = 11)
plt.ylabel(r'$\mathrm{avg.\ recall}$', fontsize = 11)
# plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.ylim([0, 1])
plt.grid()  #just add this
plt.savefig(f"{outdir}kashif_avg_recall.png")
plt.close()