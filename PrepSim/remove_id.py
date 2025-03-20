import sys

file = sys.argv[1]

assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".NoID.bed"), 'w')
for line in o:
    l = line.split("\t")
    if len(l) > 4:
        out.write("\t".join(l[0:3] + [l[4]]))
out.close()
