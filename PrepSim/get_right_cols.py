import os

for file in os.listdir():
    if "_PhyloP_Var_Dedup.bed" in file and "All" in file:
        o = open(file)
        out = open(file.replace(".bed", "_Fix.bed"), 'w')
        for line in o:
            l = line.split("\t")
            length = len(l)
            if l[length-1] == "1\n":
                out.write("\t".join(l[0:length - 6] + [l[-2]]) + "\n")
        o.close()
        out.close()