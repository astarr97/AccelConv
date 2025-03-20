### assumes sorted how bedtools requires, then by gene name ###

import sys

file = sys.argv[1]
assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".dedup.bed"), 'w')

prev = 0
for line in o:
    l = line.split("\t")
    cur = l[0] + ":" + l[2]
    if cur != prev:
        out.write(line)
    prev = cur
o.close()
out.close()