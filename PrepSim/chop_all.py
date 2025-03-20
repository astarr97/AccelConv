import sys

file = sys.argv[1]
chop_size = int(sys.argv[2])
overlap_frac = float(sys.argv[3])
overlap = overlap_frac*chop_size
assert(".bed" in file)
o = open(file)
out = open(file.replace(".bed", "_Chop" + str(chop_size) + "-" + str(overlap_frac) + ".bed"), 'w')

for line in o:
    l = line.split("\t")
    start = int(l[1])
    end = int(l[2])
    if end - start < overlap:
        out.write("\t".join(line.split("\t")[0:3]) + "\n")
    else:
        mult = ((end-start)//chop_size) + 1
        buff = (mult*chop_size - (end - start))/2
        start_new = (start - buff)//1
        end_new = (end + buff)//1
        while start_new < end_new:
            out.write("\t".join([l[0], str(max(int(start_new), 0)), str(int(start_new + chop_size))]) + "\n")
            start_new += overlap
out.close()
o.close()