import sys

file = sys.argv[1]

assert("TriPrelim_ToFilt.bed" in file)

o = open(file)
out = open(file.replace("TriPrelim_ToFilt.bed", "New_TriPrelim_ToFiltMore.bed"), 'w')
check = open("Bad_Stuff.bed", 'w')
for line in o:
    l = line.split("\t")
    if l[-2] == ".":
        out.write("\t".join(l[:-5]) + "\n")
    else:
        check.write(line)
out.close()
o.close()
check.close()
