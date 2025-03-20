import sys
import os

file = sys.argv[1]
#spec_focus is the reference species (ie positions are referenced to it)
spec_focus = sys.argv[2]
#This is the species to compare
spec_rel = sys.argv[3]
#This is the outgroup
spec_out = sys.argv[4]
spec_add_foc = sys.argv[5]
spec_add_rel = sys.argv[6]
spec_add_out = sys.argv[7]

try:
    no_phylop = sys.argv[8]
except:
    no_phylop = 0

#These are additional species to pull while we are at it
#spec_add_foc = "Hippopotamus_amphibius,Eubalaena_japonica,Eschrichtius_robustus,Balaenoptera_acutorostrata,Balaenoptera_bonaerensis,Kogia_breviceps,Platanista_gangetica,Mesoplodon_bidens,Ziphius_cavirostris,Inia_geoffrensis,Lipotes_vexillifer,Neophocaena_asiaeorientalis,Phocoena_phocoena,Delphinapterus_leucas,Monodon_monoceros,Tursiops_truncatus"
#spec_add_rel = "Tragulus_javanicus,Moschus_moschiferus,Bubalus_bubalis,Bos_taurus,Bos_indicus,Bos_mutus,Bison_bison,Beatragus_hunteri,Ammotragus_lervia,Hemitragus_hylocrius,Capra_hircus,Capra_aegagrus,Ovis_aries"
#spec_add_out = "Catagonus_wagneri"
if spec_add_rel != "NA" and spec_add_foc != "NA" and spec_add_out != "NA":
    spec_add = spec_add_foc + "," + spec_add_rel + "," + spec_add_out
elif spec_add_rel != "NA" and spec_add_foc != "NA":
    spec_add = spec_add_foc + "," + spec_add_rel
elif spec_add_rel != "NA":
    spec_add = spec_add_rel
elif spec_add_foc != "NA" and spec_add_out != "NA":
    spec_add = spec_add_foc + "," + spec_add_out
elif spec_add_foc != "NA":
    spec_add = spec_add_foc
elif spec_add_out != "NA":
    spec_add = spec_add_out
else:
    spec_add = ""

#Make folders
o = open(file)
os.mkdir(spec_focus)
os.mkdir(spec_focus + "/run1")
os.mkdir(spec_focus + "/All")

#Function to write script headers
def write_beg(out):
    out.write("#!/bin/bash\n")
    out.write("#SBATCH --time=120:00:00\n")
    out.write("#SBATCH -p hbfraser\n")
    out.write("#SBATCH --mem=8GB\n")
    out.write("#SBATCH -c 1\n\n")

#Open the first script
out = open(spec_focus + "/run1/" + "run1.sh", 'w')
write_beg(out)
base_sum = 0
c = 1

#Generic command to pull SNPs
if spec_add:
    com_snp = "/home/groups/hbfraser/Common_Software/cactus-bin-v2.6.13/bin/halSnps --minSpeciesForSnp 1 --refSequence REFSEQ --start START --length LENGTH --tsv OUT_FILE /scratch/users/astarr97/PhyloP/hg38.447way.hal " + spec_focus + " " + ",".join([spec_rel, spec_out]) + "," + spec_add + "\n"
else:
    com_snp = "/home/groups/hbfraser/Common_Software/cactus-bin-v2.6.13/bin/halSnps --minSpeciesForSnp 1 --refSequence REFSEQ --start START --length LENGTH --tsv OUT_FILE /scratch/users/astarr97/PhyloP/hg38.447way.hal " + spec_focus + " " + ",".join([spec_rel, spec_out]) + "\n"

#Command to filter SNPs
com_filt_snp = "python /scratch/users/astarr97/PhyloP/Feasibility/filter_variants.py FILE SPEC_FOCUS SPEC_REL SPEC_OUT\n"

#Commands to run PhyloP using a reference bed
com_phylop_refbed = "/home/groups/hbfraser/Common_Software/cactus-bin-v2.6.13/bin/halPhyloP /scratch/users/astarr97/PhyloP/hg38.447way.hal --refBed REFBED " + spec_focus + " /scratch/users/astarr97/PhyloP/fullTreeAnc239/fullTreeAnc239.100kb.mod test2.wig\n"

#Do PhyloP on everything
com_phylop_all = "/home/groups/hbfraser/Common_Software/cactus-bin-v2.6.13/bin/halPhyloP --refSequence REFSEQ --start START --length LENGTH /scratch/users/astarr97/PhyloP/hg38.447way.hal SPEC_FOCUS /scratch/users/astarr97/PhyloP/fullTreeAnc239/fullTreeAnc239.100kb.mod OUT_WIG\n"

#Convert to bed command
com_wig2bed = "wig2bed < OUT_WIG > OUT_BED\n"

#Command to intersect variants with PhyloP
com_intersect = "bedtools intersect -wao -a VARIANTS_BED -b PHYLOP_BED > OUT_BED\n"
out_prefix = "SNPS_"

#Function to write out one 500,000 bp block
def write_block(out, refseq, s_end, start, length, fout):
    out_wig = fout.replace(".tsv", "_PhyloP.wig")
    out_bed = out_wig.replace(".wig", ".bed")
    rel_bed = fout.replace(".tsv", "_" + spec_rel + ".bed")
    foc_bed = fout.replace(".tsv", "_" + spec_focus + ".bed")
    if length != 500000:
        out.write(com_snp.replace("REFSEQ", refseq).replace("START", str(start)).replace("LENGTH", str(length - 1)).replace("OUT_FILE", fout))
    else:
        out.write(com_snp.replace("REFSEQ", refseq).replace("START", str(start)).replace("LENGTH", str(length)).replace("OUT_FILE", fout))
    out.write(com_filt_snp.replace("FILE", fout).replace("SPEC_FOCUS", spec_focus).replace("SPEC_REL", spec_rel).replace("SPEC_OUT", spec_out))
    #out.write(com_phylop.replace("REFBED", fout.replace(".tsv", "_" + spec_focus + "_ForPhyloP.bed")))
    #out.write(com_phylop.replace("REFBED", fout.replace(".tsv", "_" + spec_rel + "_ForPhyloP.bed")))
    if not no_phylop:
        out.write(com_phylop_all.replace("SPEC_FOCUS", spec_focus).replace("REFSEQ", l[0]).replace("START", str(start)).replace("LENGTH", str(length)).replace("OUT_WIG", out_wig))
        out.write(com_wig2bed.replace("OUT_WIG", out_wig).replace("OUT_BED", out_bed))
        out.write(com_intersect.replace("VARIANTS_BED", rel_bed).replace("PHYLOP_BED", out_bed).replace("OUT_BED", rel_bed.replace(".bed", "_PhyloP_Var.bed")))
        out.write(com_intersect.replace("VARIANTS_BED", foc_bed).replace("PHYLOP_BED", out_bed).replace("OUT_BED", foc_bed.replace(".bed", "_PhyloP_Var.bed")))
    out.write("\n")

#Function to write out the commands to run at the end of a script
def write_end(out, run_num):
    out_tsv = "Run" + str(run_num) + "_" + "_ALL.tsv"
    #out.write("rm *" + "ForPhyloP.bed\n")
    if not no_phylop:
        out_foc = "Run" + str(run_num) + "_" + spec_focus + "_PhyloP_Var.bed\n"
        out_rel = "Run" + str(run_num) + "_" + spec_rel + "_PhyloP_Var.bed\n"
        out_p = "Run" + str(run_num) + "_" + spec_focus + "_PhyloP.bed\n"
        out.write("cat *" + spec_focus + "*PhyloP_Var.bed > ../All/" + out_foc)
        out.write("cat *" + spec_rel + "*PhyloP_Var.bed > ../All/" + out_rel)
        out.write("cat *PhyloP.bed > ../All/" + out_p)
        out.write("cat *.tsv > ../All/" + out_tsv)
    else:
        out_foc = "Run" + str(run_num) + "_" + spec_focus + "_ForPhyloP_Var.bed\n"
        out_rel = "Run" + str(run_num) + "_" + spec_rel + "_ForPhyloP_Var.bed\n"
        out.write("cat *" + spec_focus + "*ForPhyloP_Var.bed > ../All/" + out_foc)
        out.write("cat *" + spec_rel + "*ForPhyloP_Var.bed > ../All/" + out_rel)
        out.write("cat *.tsv > ../All/" + out_tsv)
    out.write("\n")

run_contigs = 0
for line in o:
    run_contigs += 1
    l = line.replace("\n", "").split("\t")
    cur_len = int(l[2])
    #Iterate through in batches of 500,000 bases, resetting and opening a new script if we reach 20,000,000 bases
    if cur_len < 500000:
        fout = out_prefix + l[0] + "_" + str(0) + "-" + str(int(l[2])) + ".tsv"
        write_block(out, l[0], int(l[2]), 0, int(l[2]), fout)
        base_sum += cur_len
    else:
        s_end = 500000
        while s_end < cur_len:
            if base_sum > 20000000:
                write_end(out, c)
                out.close()
                c += 1
                os.mkdir(spec_focus + "/" + "run" + str(c))
                out = open(spec_focus + "/" + "run" + str(c) + "/" + "run" + str(c) + ".sh", 'w')
                write_beg(out)
                base_sum = 0
                run_contigs = 0
            fout = out_prefix + l[0] + "_" + str(s_end - 500000) + "-" + str(s_end) + ".tsv"
            write_block(out, l[0], s_end, s_end - 500000, 500000, fout)
            base_sum += 500000
            s_end += 500000
        fout = out_prefix + l[0] + "_" + str(s_end - 500000) + "-" + str(l[2]) + ".tsv"
        write_block(out, l[0], s_end, s_end - 500000, int(l[2]) - (s_end - 500000), fout)
        base_sum += int(l[2]) - (s_end - 500000)
        

    if base_sum > 20000000 or run_contigs >= 300:
        write_end(out, c)
        out.close()
        c += 1
        os.mkdir(spec_focus + "/" + "run" + str(c))
        out = open(spec_focus + "/" + "run" + str(c) + "/" + "run" + str(c) + ".sh", 'w')
        write_beg(out)
        base_sum = 0
        run_contigs = 0
write_end(out, c)
out.close()
