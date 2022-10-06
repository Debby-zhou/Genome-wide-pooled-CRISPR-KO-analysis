import pandas as pd
import numpy as np
from sys import argv

counts_file = "sgrna_counts/geckov2_sgrnas.count_normalized.txt"
rivals = "sample1_7day,sample1_21day"
output_file = "gene-based_log2FC_sample1.txt"

arrRivals = [r.split(",") for r in rivals.split(";")] if ";" in rivals else [rivals.split(",")]

def log2fc(early_time, late_time):
    log2fcNum = np.log2(late_time+1)-np.log2(early_time+1)
    return log2fcNum

data = pd.read_csv(counts_file, sep="\t")
times = 0
for team in arrRivals:
	early, late = team[0], team[1]
	df_rival = data[['sgRNA','Gene']+team]
	log2fcName = "log2FC("+late+"/"+early+")"
	absLog2fcName = "abs_"+log2fcName
	df_rival[log2fcName] = df_rival.apply(lambda x: log2fc(x[early],x[late]), axis=1)
	df_rival[absLog2fcName] = abs(df_rival[log2fcName])
	df_result = df_rival.sort_values(by=absLog2fcName, ascending=False).drop_duplicates(subset='Gene')
	if times==0:
		df_final = df_result
	else:
		df_final = pd.merge(df_final,df_result,how='inner',left_on="Gene",right_on="Gene")
	times += 1
df_final.to_csv(output_file, sep="\t", index=False)
