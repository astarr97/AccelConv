import sys

file = sys.argv[1]

assert(".fasta" in file)

o = open(file)
out = open(file.replace(".fasta", ".fastabed"), 'w')

to_write = []
for line in o:
    if ">" in line:
        l = line.replace("\n", "").replace(">", "")
        to_write = [l.split(":")[0], str(int(l.split(":")[1].split("-")[0]) + 1), str(int(l.split(":")[1].split("-")[1]) - 1)]
    else:
        to_write.append(line)
        out.write("\t".join(to_write))
o.close()
out.close()
