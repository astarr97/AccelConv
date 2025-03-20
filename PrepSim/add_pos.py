import sys

file = sys.argv[1]

assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".pos.bed"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join(l + [l[0] + ":" + l[2]]) + "\n")
o.close()
out.close()
