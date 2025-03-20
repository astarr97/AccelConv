import sys

file = sys.argv[1]
phylo_col = int(sys.argv[2])
gene_col = int(sys.argv[3])
old_tri_col = int(sys.argv[4])
try:
    stay_bed = int(sys.argv[5])
except:
    stay_bed = 0

assert("AncPrelim.bed" in file)

o = open(file)
out = open(file.replace("Prelim", ""), 'w')

for line in o:
    l = line.split("\t")
    if l[-2] != ".":
        if stay_bed:
            out.write("\t".join([l[0], l[1], l[2]] + [l[phylo_col], l[gene_col], l[-2].upper()]) + "\n")
        else:
            out.write("\t".join([l[0] + ":" + l[2]] + [l[phylo_col], l[gene_col], l[-2].upper()]) + "\n")
    else:
        if stay_bed:
            out.write("\t".join([l[0], l[1], l[2]] + [l[phylo_col], l[gene_col], l[old_tri_col].upper()]) + "\n")
        else:
            out.write("\t".join([l[0] + ":" + l[2]] + [l[phylo_col], l[gene_col], l[old_tri_col].upper()]) + "\n")
o.close()
out.close()
