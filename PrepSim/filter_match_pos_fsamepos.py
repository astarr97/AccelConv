import sys

file = sys.argv[1]
out_file = sys.argv[2]

o = open(file)
out = open(out_file, 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    if (l[0] == l[-5].split(":")[0] and int(l[2]) == int(l[-5].split(":")[1])) and l[-1] == "0":
        out.write("\t".join(l[0:-5]) + "\n")
o.close()
out.close()
