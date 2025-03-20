import sys

gtf = sys.argv[1]
name_keep = sys.argv[2]

assert(".gtf" in gtf)
o = open(gtf)
out = open(gtf.replace(".gtf", "_CDS.bed"), 'w')
out2 = open(gtf.replace(".gtf", "_Start.bed"), 'w')
for line in o:
    if "#" not in line:
        l = line.split("\t")
        if l[2] == "CDS":
            ll = l[8].split(";")
            for i in ll:
                if name_keep in i:
                    gene = i.replace('"', '').split(" ")[1].split(".")[1]
            out.write("\t".join([l[0], str(int(l[3]) - 1), l[4], gene]) + "\n")
        #elif l[2] == "start_codon":
        #    ll = l[8].split(";")
        #    for i in ll:
        #        if name_keep in i:
        #            gene = i.replace('"', '').split(" ")[1].split(".")[1]
        #    out2.write("\t".join([l[0], str(int(l[3]) - 1), l[4], gene]) + "\n")
out.close()
o.close()
