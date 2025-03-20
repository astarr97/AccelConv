import sys

file = sys.argv[1]
pos = int(sys.argv[2])

assert(".bed" in file)

o = open(file)
out_cds = open(file.replace(".bed", ".CDS.bed"), 'w')
out_nc = open(file.replace(".bed", ".NonCod.bed"), 'w')
for line in o:
    l = line.replace("\n", "").split("\t")
    if int(l[len(l)-1]) == 0:
        out_cds.write("\t".join(l[0:pos] + l[len(l)-2:len(l)]) + "\n")
    elif int(l[len(l)-1]) > 0:
        out_nc.write("\t".join(l[0:pos] + l[len(l)-2:len(l)]) + "\n")
out_cds.close()
out_nc.close()
