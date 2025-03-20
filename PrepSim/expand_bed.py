import sys

file = sys.argv[1]
chrom_size = sys.argv[2]

d = {}
o = open(chrom_size)
for line in o:
    l = line.replace("\n", "").split("\t")
    d[l[0]] = int(l[2])

assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".exp.bed"), 'w')

for line in o:
    if "#" != line[0]:
        l = line.split("\t")
        if (int(l[1]) - 1) >= 0 and (int(l[2]) + 1) <= d[l[0]]:
            out.write("\t".join([l[0], str(int(l[1]) - 1), str(int(l[2]) + 1)]) + "\n")
o.close()
out.close()
