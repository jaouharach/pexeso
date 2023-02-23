from os import listdir
from os.path import isfile, join
import pandas as pd


k_values = ["1", "5", "10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000"]
k_values = ["1000000"]
src_dir_path = "/home/jaouhara.chanchaf/work-dir/exp-results/kashif-search/100tk/stop-mode-0/kashif_l100000_10q_min100_max100/"


def merge_files(dir_path, k): # merge result file for a given k
    result_files = [f for f in listdir(dir_path) if isfile(join(dir_path, f))]
    result_list = []
    
    for i, f in enumerate(result_files):
        if f"_k{k}_" in f:
            result_list.append(pd.read_csv(f"{dir_path}/{f}"))
    
    print(result_list)
    results_merged = pd.concat(result_list, ignore_index=True)
    results_merged.to_csv(f'./merged_files/merged_results_k={k}.csv', index=False)


for k in k_values:
    merge_files(src_dir_path, k)
