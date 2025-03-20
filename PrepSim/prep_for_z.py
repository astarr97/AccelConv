import sys
import pandas as pd
import numpy as np
import os

folder_name = sys.argv[1]
prefix = sys.argv[2]

folder = os.listdir(folder_name)
folder.sort(key = lambda x: int(x.split("_")[0].replace("Simulation", "")))

#start_cols = ["Species1 Sum PhyloP", "Num Species1 Var", "Species2 Sum PhyloP", "Num Species2 Var", "PhyloP Difference", "PhyloP L2FC", "Median_PhyloP_Species1", "Median_PhyloP_Species2", "MWU p-value", "Centered PhyloP L2FC", "Centered PhyloP Difference"]
#start_cols = ["Species1 Sum PhyloP", "Num Species1 Var", "Species2 Sum PhyloP", "Num Species2 Var", "PhyloP Difference", "PhyloP L2FC", "Centered PhyloP L2FC", "Centered PhyloP Difference"]
start_cols = ["Species1 Sum PhyloP", "Num Species1 Var", "Species2 Sum PhyloP", "Num Species2 Var", "PhyloP Difference", "PhyloP L2FC"]

add_cols = ['PhyloP Difference Sim', 'PhyloP L2FC Sim', 'Num Species1 Var Sim', 'Num Species2 Var Sim']
ind = 1
v_start = 0
v_temp = 0
seen = []
c = 0
for file in folder:
    c += 1
    if c % 1000 == 0:
        print(c)
    v = pd.read_csv(folder_name + "/" + file, sep = "\t")
    v = v.set_index(v.columns[0])
    if ind:
        v_start = v[start_cols].copy()
        v_temp = v[[]]
        ind = 0
    v = v[add_cols].copy()
    for col in add_cols:
        if "Var" not in col:
            v[col] = np.round(v[col], 3)
    seen.append(file.split("_")[0])
    v.columns = [x + " " + file.split("_")[0] for x in add_cols]
    v_temp = v_temp.join(v)

v_temp = v_temp.astype(str)
for col in add_cols:
    out_col = []
    v_col = v_temp[[col + " " + x for x in seen]]
    for index, row in v_col.iterrows():
        out_col.append([index, ";".join(list(row))])
    to_join = pd.DataFrame(out_col).set_index(0)
    to_join.columns = [col]
    v_start = v_start.join(to_join)
v_start.to_csv(prefix + "AllSims.txt", sep = "\t")
