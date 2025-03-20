import sys

file = sys.argv[1]

o = open(file)
prev_gene = 0
for line in o:
    l = line.split("\t")
    gene = l[3].replace("\n", "")
    if prev_gene != gene:
        if prev_gene:
            out.close()
        out = open("Dedup_CDS/" + gene + "_CDS.bed", 'w')
    out.write(line)
    prev_gene = gene
o.close()
out.close()
