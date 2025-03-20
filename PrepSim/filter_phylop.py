import sys


file = sys.argv[1]
cutoff = int(sys.argv[2])
pos = int(sys.argv[3])

assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", "_FiltCut" + str(cutoff) + ".bed"), 'w')

for line in o:
    l = line.split("\t")
    if float(l[pos]) >= cutoff:
        out.write(line)
o.close()
out.close()
