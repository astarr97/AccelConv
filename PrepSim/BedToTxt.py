import sys
import numpy as np

file = sys.argv[1]
assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".txt"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join([l[0] + ":" + l[2]] + l[3:]) + "\n")
o.close()
out.close()
