import sys
import pandas as pd

chrom_file = sys.argv[1]
bed = sys.argv[2]

assert(".bed" in bed)

v = pd.read_csv(chrom_file, sep = "\t")

d_conv = {}
for index, row in v.iterrows():
    d_conv[row["UCSC style name"]] = row["GenBank seq accession"]

o = open(bed)
out = open(bed.replace(".bed", ".ConvChrom.bed"), 'w')

for line in o:
    l = line.split("\t")
    try:
        l[0] = d_conv[l[0]]
        out.write("\t".join(l))
    except:
        print(l)
o.close()
out.close()
