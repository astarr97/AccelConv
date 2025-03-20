import sys
import pandas as pd

to_filt = sys.argv[1]
bad = sys.argv[2]
gene_col = int(sys.argv[3])

assert("ToFiltMore.bed" in to_filt)

v = pd.read_csv(to_filt, sep = "\t", header = None)
blacklist = list(pd.read_csv(bad, sep = "\t", header = None)[0])
v = v[~v[gene_col].isin(blacklist)]
v.to_csv(to_filt.replace("_ToFiltMore", ""), sep = "\t", header = False, index = False)
