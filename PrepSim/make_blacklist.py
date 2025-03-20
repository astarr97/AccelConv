import sys
from collections import Counter
import os

to_filt = sys.argv[1]
bad = sys.argv[2]
gene_col = int(sys.argv[3])
filt_prop = float(sys.argv[4])
blacklist_file = sys.argv[5]

assert("ToFiltMore.bed" in to_filt)


o = open(bad)
#out = open(to_filt.replace("_ToFiltMore.bed", ".bed"), 'w')

bad_counter = {}
for line in o:
    l = line.split("\t")
    bad_counter[l[gene_col]] = 0
o.close()

o = open(bad)
for line in o:
    l = line.split("\t")
    bad_counter[l[gene_col]] += 1
o.close()

o = open(to_filt)
pass_counter = {}
for line in o:
    l = line.split("\t")
    pass_counter[l[gene_col]] = 0
o.close()

o = open(to_filt)
for line in o:
    l = line.split("\t")
    pass_counter[l[gene_col]] += 1
o.close()

if blacklist_file in os.listdir():
    badder = open(blacklist_file, 'a')
else:
    badder = open(blacklist_file, 'w')
for key in bad_counter.keys():
    if key not in pass_counter.keys():
        badder.write(key + "\n")
    #We divide by 2 as all sites are double counted due to overlapping windows
    elif bad_counter[key]/(bad_counter[key] + pass_counter[key])/2 >= filt_prop:
        badder.write(key + "\n")

#Remove Cadherins
for key in pass_counter.keys():
    if "PCDHG" in key or "PCDHA" in key or "PCDHB" in key or "UGT1A" in key:
        badder.write(key + "\n")
badder.close()
