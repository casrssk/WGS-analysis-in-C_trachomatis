## importing libraries
import sys
import os
import subprocess
import pandas as pd
import numpy as np
dir_in = "results/var_call/filtered/"
dir_vcf = "results/var_call/"


def vcf_com(vcf1, vcf2, snp_table1, snp_table2):
	
	filepath1 = os.path.join(dir_vcf, sample1)
	filepath2 = os.path.join(dir_vcf, sample2)
	filepath3 = os.path.join(dir_in, snp_table1)
	filepath4 = os.path.join(dir_in, snp_table2)
	snp_concordance = "results/var_call/snp_concordance.txt"
	compare = "  ".join(["SnpSift concordance -v", filepath1, filepath2, ">", snp_concordance])
	subprocess.call(compare, shell=True)
	snp1 = pd.read_csv(filepath3, sep='\t')
	snp2 = pd.read_csv(filepath4, sep='\t')
	annot_dict1 = dict(zip(snp1["POS"], snp1["ANN"]))
	annot_dict2 = dict(zip(snp2["POS"], snp2["ANN"]))
	snp_vis = pd.merge(snp1, snp2, how='outer', on=["POS"])
	snp_vis = snp_vis.sort_values("POS")
	snp_vis = snp_vis[["POS", "GT_x", "GT_y"]]
	snp_vis = snp_vis.fillna(0)
	snp_vis = snp_vis.replace(to_replace=r'^[AGCT]$', value=1, regex=True)
	diff = snp_vis.loc[(snp_vis["GT_x"] == 0) | (snp_vis["GT_y"] == 0)]
	diff = diff.reset_index(drop=True)
	print(len(diff))
    snp_vis = snp_vis.set_index('POS')
    fig, ax = plt.subplots(figsize=(8, 18))
    sns_plot = sns.heatmap(snp_vis, ax=ax)
    figure = sns_plot.get_figure()
    figure.savefig("results/var_call/vcf_compare.png", dpi=400)
    for i in range(len(diff)):
        if diff["GT_x"][i] == 1:
            for k, v in annot_dict1.items():
                if k == diff['POS'][i]:
                    diff['annot'] = v
        else:
            for k, v in annot_dict2.items():
                if k == diff['POS'][i]:
                    diff["annot"] = v
    return diff
