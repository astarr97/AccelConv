import sys

file = sys.argv[1]

assert(".bed") in file
o = open(file)
out = open(file.replace(".bed", ".point1.bed"), 'w')

for line in o:
    l = line.split("\t")
    l[0] = l[0] + ".1"
    out.write("\t".join(l))
o.close()
out.close()
