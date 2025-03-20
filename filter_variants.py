import pandas as pd
import sys
import numpy as np

file = sys.argv[1]
spec_focus = sys.argv[2]
spec_rel = sys.argv[3]
spec_out = sys.argv[4]

assert(".tsv" in file)

#Read in the file
v = pd.read_csv(file, sep = "\t")
v = v.replace(np.nan, "N")
#Convert all variants to upper case
for col in list(v.columns):
    if "refSequence" != col and "refPosition" != col:
        v[col] = v[col].str.upper()
        
#Write the file out again, replacing it now that everything is uppercase
v.to_csv(file, sep = "\t", index = False)

#Filter out sites that are N or the same in spec_focus and spec_rel
v = v[(v[spec_focus] != v[spec_rel]) & (v[spec_focus].isin(["A", "C", "G", "T"])) & (v[spec_rel].isin(["A", "C", "G", "T"])) & (v[spec_out].isin(["A", "C", "G", "T"]))]

#Convert to bed format
v["refPosition1"] = v["refPosition"] + 1
v = v[["refSequence", "refPosition", "refPosition1", spec_focus, spec_rel, spec_out]]
v = v.dropna()

#Split into spec_focus-derived and spec_rel-derived
v_focus_der = v[v[spec_rel] == v[spec_out]]
v_rel_der = v[v[spec_focus] == v[spec_out]]

#Write out to bed files
v_focus_der.to_csv(file.replace(".tsv", "_" + spec_focus + ".bed"), sep = "\t", header = False, index = False)
v_rel_der.to_csv(file.replace(".tsv", "_" + spec_rel + ".bed"), sep = "\t", header = False, index = False)

#Unused 
#v_focus_der[["refSequence", "refPosition", "refPosition1"]].to_csv(file.replace(".tsv", "_" + spec_focus + "_ForPhyloP.bed"), sep = "\t", header = False, index = False)
#v_rel_der[["refSequence", "refPosition", "refPosition1"]].to_csv(file.replace(".tsv", "_" + spec_rel + "_ForPhyloP.bed"), sep = "\t", header = False, index = False)
