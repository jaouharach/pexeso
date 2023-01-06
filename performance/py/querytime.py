import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

base_path = Path(__file__).parent
outdir = '../img/'
csv_files = {"kashif": (base_path / "../csv/kashif_querytime.csv").resolve(), "pexeso": (base_path / "../csv/pexeso_querytime.csv").resolve()}

# plt.rcParams["figure.figsize"] = (8,7)
dfk = pd.read_csv(csv_files['kashif'])
dfp = pd.read_csv(csv_files['pexeso'])

print(dfk)
print(dfp)
