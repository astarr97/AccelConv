from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import numpy as np
import sys

sims_file = sys.argv[1]
out_file = sys.argv[2]
sims = pd.read_csv(sims_file, sep = "\t")
sims = sims.set_index(sims.columns[0])

#exclude = []
#for index, row in sims.iterrows():
#    if "PCDHG" in index and "PCDHGA2" not in index:
#        exclude.append(index)
#sims = sims.loc[np.setdiff1d(sims.index, exclude)]

zscore_dif = []
zscore_l2fc = []
pvalue_dif = []
pvalue_l2fc = []
for index, row in sims.iterrows():
    difs = [np.float64(x) for x in row['PhyloP Difference Sim'].split(";")]
    l2fcs = [np.float64(x) for x in row['PhyloP L2FC Sim'].split(";")]
    stdev_dif = np.std(difs)
    mean_dif = np.mean(difs)
    stdev_l2fc = np.std(l2fcs)
    mean_l2fc = np.mean(l2fcs)
    z_dif = (row["PhyloP Difference"] - mean_dif)/stdev_dif
    z_l2fc = (row["PhyloP L2FC"] - mean_l2fc)/stdev_l2fc
    zscore_dif.append(z_dif)
    zscore_l2fc.append(z_l2fc)
    p_dif = norm.sf(abs(z_dif))*2
    p_l2fc = norm.sf(abs(z_l2fc))*2
    pvalue_dif.append(p_dif)
    pvalue_l2fc.append(p_l2fc)

sims_write = sims[["Species1 Sum PhyloP", "Num Species1 Var", "Species2 Sum PhyloP", "Num Species2 Var", "PhyloP Difference", "PhyloP L2FC"]].copy()

sims_write["Z-score Difference"] = zscore_dif
sims_write["Z-score L2FC"] = zscore_l2fc
sims_write["p-value Difference"] = pvalue_dif
sims_write["p-value L2FC"] = pvalue_l2fc
sims_write["FDR Difference"] = fdrcorrection(pvalue_dif)[1]
sims_write["FDR L2FC"] = fdrcorrection(pvalue_l2fc)[1]
sims_write = sims_write.sort_values("FDR L2FC")
sims_write.to_csv(out_file, sep = "\t")
