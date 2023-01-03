import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (6,5)
df = pd.read_csv('stats.csv')
# print(df)
df['tau'] = 100. * df['tau']

dfc = df[['%filled_cells', 'tau', '#visited_cells', '#filtered_cells', '#visited_mcells', '#visited_ccells']].copy()
dfv = df[['%filled_cells', 'tau', '#filtered_vectors', '#checked_pvectors', '#checked_mvectors']].copy()

# box plot to show impact of data distribution in the pivot space on the filtering power of pexeso
dfvv = dfv.melt(['tau', '%filled_cells'], var_name='metric', value_name='val')
print(dfvv)
g = sns.boxplot(data=dfvv, x='%filled_cells', y='val', hue='metric', palette=sns.color_palette('icefire', 10))
plt.xlabel(r'$\mathrm{\%\ of\ non\ empty\ leaf\ cells}$', fontsize = 11)
plt.ylabel(r'$\mathrm{value}$', fontsize = 11)
plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"impact_of_data_distribution_on_data_access.png")
plt.close()


# box plot to show impact of data distribution in the pivot space on accessing cells in pexeso grid
dfcc = dfc.melt(['tau', '%filled_cells'], var_name='metric', value_name='val')
# print(dfcc)
g = sns.boxplot(data=dfcc, x='%filled_cells', y='val', hue='metric', palette=sns.color_palette("rocket", 4))
plt.xlabel(r'$\mathrm{\%\ of\ non\ empty\ leaf\ cells}$', fontsize = 11)
plt.ylabel(r'$\mathrm{value}$', fontsize = 11)
plt.legend(loc='best')
# plt.yscale("log")
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"impact_of_data_distribution_on_cell_access.png")
plt.close()
