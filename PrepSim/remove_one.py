import sys

file = sys.argv[1]

assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".NoOne.bed"), 'w')
for line in o:
    l = line.split("\t")
    out.write("\t".join([l[0]] + [l[2]]))
out.close()
