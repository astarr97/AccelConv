import sys

file = sys.argv[1]
out_file = sys.argv[2]

o = open(file)
out = open(out_file, 'w')

prev = 0
for line in o:
    l = line.replace("\n", "").split("\t")
    if l[-1] == prev:
        out.write("\t".join([l[-1].split(":")[0], str(int(l[-1].split(":")[1])-1), l[-1].split(":")[1]]) + "\n")
    prev = l[-1]
o.close()
out.close()
