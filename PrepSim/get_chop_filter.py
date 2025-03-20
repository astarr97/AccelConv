import sys

file = sys.argv[1]
cutoff = int(sys.argv[2])

assert("_Count.bed" in file)

o = open(file)
out = open(file.replace("_Count.bed", "_ToSubtract.bed"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    if int(l[3]) >= cutoff:
        out.write(line)
o.close()
out.close()
