import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

outdir = '../img/'
csv_file = '../nn_insert_time.csv'

plt.rcParams["figure.figsize"] = (8,7)
df = pd.read_csv(csv_file)

df['k'] = df['k'].astype('str')
g = sns.lineplot(data=df, x='k', y='total_insert_time')

# g.set_xticklabels(rotation=30)
plt.subplots_adjust(bottom=0.15)
plt.xlabel(r'$\mathrm{k}$', fontsize = 14)
plt.ylabel(r'$\mathrm{total\ nn\ insert\ time}$', fontsize = 14)
# plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"{outdir}total_insert_time.png")
plt.close()



