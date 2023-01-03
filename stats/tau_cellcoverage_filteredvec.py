import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

outdir = './img-tau[0-100]/'
csv_file = 'stats.csv'

plt.rcParams["figure.figsize"] = (6,5)
df = pd.read_csv(csv_file)
# print(df)
df['tau'] = 100. * df['tau']


# plot change in #filtered_vector for different %filled_cells
dfv = df[['%filled_cells', 'tau', '#filtered_vectors']].copy()
dfvv = dfv.melt(['tau', '%filled_cells'], var_name='metric', value_name='val')
print(dfvv)

g = sns.catplot(data=dfvv, x='tau', y='val', hue='%filled_cells',  kind='point', scale = 0.5,
    palette=sns.color_palette('icefire', 4),  markers=['o', 'v', '*', 'X'])
g.set_xticklabels(rotation=30)
plt.xlabel(r'$\mathrm{tau}$', fontsize = 11)
plt.ylabel(r'$\mathrm{nb\ filtered\ vectors}$', fontsize = 11)
plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"{outdir}cellcoverage_tau_on_filtered_vectors.png")
plt.close()


# plot change in #filtered_cells for different %filled_cells
dfc = df[['%filled_cells', 'tau', '#filtered_cells']].copy()
dfcc = dfc.melt(['tau', '%filled_cells'], var_name='metric', value_name='val')
print(dfcc)

g = sns.catplot(data=dfcc, x='tau', y='val', hue='%filled_cells',  kind='point', scale = 0.5,
    palette=sns.color_palette('icefire', 4),  markers=['o', 'v', '*', 'X'])
g.set_xticklabels(rotation=30)
plt.xlabel(r'$\mathrm{tau}$', fontsize = 11)
plt.ylabel(r'$\mathrm{nb\ filtered\ cells}$', fontsize = 11)
plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"{outdir}cellcoverage_tau_on_filtered_cells.png")
plt.close()

