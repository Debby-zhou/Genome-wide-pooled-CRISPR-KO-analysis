import numpy as np
import pandas as pd
from sys import argv
'''
nor_counts_file = "sgrna_counts/geckov2_sgrnas.count_normalized.txt"
fc_samples = "sample1_7day,sample1_21day" # with multiple samples please use ';' to seperate e.g. "sample1_7day,sample1_21day;sample2_7day,sample2_21day"
output_file = "gene-based_log2FC_sample1.txt"
'''
script, nor_counts_file, fc_samples, output_file = argv 
fc_rivals = [r.split(",") for r in fc_samples.split(";")] if ";" in fc_samples else [fc_samples.split(",")]
data = open(nor_counts_file).readlines()
header = data[0].rstrip("\n").split("\t") # 0: sgRNA, 1: gene, 2~: samples
data = [line.rstrip("\n").split("\t") for line in data[1:]]
arr_expr = [[float(split_line[split_index]) for split_index in range(len(split_line)) if split_index>1] for split_line in data]
arr_expr = np.array(arr_expr)
df_final = pd.DataFrame()
for rival in fc_rivals:
	early_sample, late_sample = rival[0], rival[1]
	early_key = [h for h in range(len(header)) if header[h]==early_sample][0]-2
	late_key = [h for h in range(len(header)) if header[h]==late_sample][0]-2
	log2fc =[np.log2(row[late_key]+1)-np.log2(row[early_key]+1) for row in arr_expr]
	log2fc_name = "log2("+late_sample+"/"+early_sample+")"
	df_log2fc = pd.DataFrame({
		header[1]: [d[1] for d in data],
		header[0]+log2fc_name[4:]: [d[0] for d in data],
		early_sample: [e[early_key] for e in arr_expr],
		late_sample: [e[late_key] for e in arr_expr],
		log2fc_name: log2fc,
		"abs_"+log2fc_name: np.absolute(log2fc)
	})
	df_result = df_log2fc.sort_values(by="abs_"+log2fc_name, ascending=False).drop_duplicates(subset=header[1])
	if df_final.empty:
		df_final = df_result
	else:
		df_final = pd.merge(df_final,df_result,how='inner',left_on=header[1],right_on=header[1]) 
df_final.to_csv(output_file, index=False) if ".csv" in output_file else df_final.to_csv(output_file, sep="\t", index=False)
