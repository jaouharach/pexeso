# pexeso results are all stored in one file
# kashif results are stored in several files each file stores the results of one query column (for a given k values)

from os import listdir
from os.path import isfile, join
import pandas as pd
import csv 


def format_query_id(query_id): # from [tableid]:[columnid] to t[tableid]c[columnid]
    tableid, setid = query_id.strip().split(":") 
    return f"t{tableid}c{setid}"

# k_values = ["1", "5", "10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000"]
k_values = ["1000000"]

kashif_results_path = "./merged_files/"
pexeso_file = "./100k_tables_10queries_size100/queryresults_tau=6%_T=60%.csv"

out_file = './query_recall_k=1m.csv'
f = open(out_file, 'w')
writer = csv.writer(f)
writer.writerow(["tqq","recall","precision", "k"])

pexeso_df = pd.read_csv(pexeso_file)
pexeso_df = pexeso_df.drop('nooverlap', axis=1)

for k in k_values:
    
    kashif_file = [kf for kf in listdir(kashif_results_path) if isfile(join(kashif_results_path, kf)) and f'_k={k}.' in kf][0]
    
    kashif_df = pd.read_csv(f"{kashif_results_path}{kashif_file}")
    kashif_df = kashif_df.drop([' q', ' s', ' q_pos', ' s_pos', ' time', ' k', ' d'], axis=1)

    queryids = pd.unique(kashif_df['TQ:Q'])

    for query in queryids:
        print(f"processing query: ({format_query_id(query)})")

        res = kashif_df.loc[kashif_df['TQ:Q'] == query].reset_index(drop=True) # kashif results
        gt = pexeso_df.loc[pexeso_df['tqq'] == format_query_id(query)].reset_index(drop=True) # ground truth results (pexeso)

        res['TQ:Q'] = res['TQ:Q'].apply(lambda qid : format_query_id(qid))
        res[' TS:S'] = res[' TS:S'].apply(lambda qid : format_query_id(qid))

        res['tp'] = res[' TS:S'].map(gt['tss'].value_counts())
        res['tp'] = res['tp'].fillna(0)
        

        tmpr = pd.unique(res[' TS:S'])
        tmpgt = pd.unique(gt['tss'])
        tp = list(set(tmpr) & set(tmpgt))

        # print("results")
        # print(tmpr)
        # print("ground truths")
        # print(tmpgt)
        # print("TPs")
        # print(tp)

        if gt.shape[0] != 0:
            recall = len(tp)/ len(tmpgt)
            precision = len(tp)/ len(tmpr)
            writer.writerow([format_query_id(query), recall, precision, k])
            print(f"query ({query}), k = {k}, recall = {recall}, precision = {precision}")


    
