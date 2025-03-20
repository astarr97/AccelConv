import sys

file = sys.argv[1]
phylo_col = int(sys.argv[2])
gene_col = int(sys.argv[3])
is_focal = int(sys.argv[4])
try:
    stay_bed = int(sys.argv[5])
except:
    stay_bed = 0

assert("TriPrelim.bed" in file or "interback.bed" in file)

o = open(file)
out = open(file.replace("Prelim", "").replace("interback", "toinp"), 'w')

d_comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
prev = 0
for line in o:
    l = line.split("\t")
    if l[-2] != "." and l[0] + ":" + l[2] != prev:
        if is_focal:
            new_center = l[4]
            if l[-2][1].upper() != l[3].upper():
                new_center = d_comp[new_center.upper()]
            ancestral = l[-2][0] + new_center + l[-2][2]
            if stay_bed:
                out.write("\t".join([l[0], l[1], l[2]] + [l[phylo_col], l[gene_col], l[-2], ancestral]) + "\n")
            else:
                out.write("\t".join([l[0] + ":" + l[2]] + [l[phylo_col], l[gene_col], l[-2], ancestral]) + "\n")
        else:
            if stay_bed:
                out.write("\t".join([l[0], l[1], l[2]] + [l[phylo_col], l[gene_col], l[-2]]) + "\n")
            else:
                out.write("\t".join([l[0] + ":" + l[2]] + [l[phylo_col], l[gene_col], l[-2]]) + "\n")
    prev = l[0] + ":" + l[2]
o.close()
out.close()
