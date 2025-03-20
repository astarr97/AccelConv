import sys
import pandas as pd
import numpy as np
import os
from collections import Counter
from copy import deepcopy

species1_file = sys.argv[1]
species2_file = sys.argv[2]
background_file = sys.argv[3]

### NOTE ###
#It does not matter that we could never have a central A or T trinucleotide and PhyloP > 7.
#It will never come up later so we can assign a probability value with no issues.

#Maximum allowable PhyloP score for non-background
#Typically 8.5 for our data
#Higher values indicate misalignment and/or other funny business
#Even though A/T nucleotides cannot generally have PhyloP scores this big, A/T nucleotides with PhyloP scores greater than 7.43 (the equivalent cutoff) do not indicate problems
#This is because they more likely indicate a nucleotide that was ancestrally A/T but became G/C in the common ancestor of the two species
max_phylo = float(sys.argv[4])

#Folder to write out the tables to
folder = sys.argv[5]

#What to control for during the simulation when assigning variants to species and gene
#Can be one of three options:
#(1) just_tri means only control for the trinucleotide context
#(2) tri_species_score means additionally controlling for the score (often phylop, but not necessarily) only when assigning to species
#(3) tri_both_score means additionally controlling for the score when assigning to both species and gene (not possible unless score is computed for all background variants)
#(4) just_mono means only control for the mononucleotide context (typically used for more closely related species)
#(5) mono_species_score is analogous to (2)
#(6) mono_both_score is analogous to (3)
control = sys.argv[6]

#Minimum number of variants needed to be included
#Has format integer,and or integer,or; and implies both species must have at least integer variants; or implies just one species does
int_op = sys.argv[7]
num_var_cut = int(int_op.split(",")[0])
operation = int_op.split(",")[1]

if folder not in os.listdir():
    os.mkdir(folder)

assert("Anc.bed" in background_file and ".bed" in species1_file and ".bed" in species2_file)

phylop_col = 1
gene_col = 2
tri_col = 3
anc_tri_col = 4

species1 = pd.read_csv(species1_file, sep = "\t", header = None).drop_duplicates()
species2 = pd.read_csv(species2_file, sep = "\t", header = None).drop_duplicates()

alphabet = ["A", "T", "C", "G"]
toss = ["NNN"]
for i in alphabet:
    toss.append(i + "NN")
    toss.append("N" + i + "N")
    toss.append("NN" + i)
    for j in alphabet:
        toss.append(i + j + "N")
        toss.append(i + "N" + j)
        toss.append("N" + i + j)
        toss.append(j + i + "N")
        toss.append(j + "N" + i)
        toss.append("N" + j + i)
toss = list(set(toss))

species1[tri_col] = species1[tri_col].str.upper()
species2[tri_col] = species2[tri_col].str.upper()
species1[anc_tri_col] = species1[anc_tri_col].str.upper()

if "mono" in control:
    species1[tri_col] = species1[tri_col].str[1]
    species2[tri_col] = species2[tri_col].str[1]
    species1[anc_tri_col] = species1[anc_tri_col].str[1]
    toss = ["N"]

#Some preliminary filtering
species1 = species1[~species1[phylop_col].isin(["."])]
species2 = species2[~species2[phylop_col].isin(["."])]
species1 = species1[~species1[tri_col].isin(toss)]
species2 = species2[~species2[tri_col].isin(toss)]
species1[phylop_col] = species1[phylop_col].astype(float)
species2[phylop_col] = species2[phylop_col].astype(float)

#Usually superfluous, but we don't want negative PhyloP to count against genes
species1[phylop_col] = np.maximum(species1[phylop_col], 0)
species2[phylop_col] = np.maximum(species2[phylop_col], 0)

#Remove PhyloP scores that we know are wrong
species1 = species1[species1[phylop_col] <= max_phylo]
species2 = species2[species2[phylop_col] <= max_phylo]

print(len(np.unique(list(species1[gene_col]) + list(species2[gene_col]))))

#Get list of all genes that will be included
num_var_s1 = Counter(species1[gene_col])
num_var_s2 = Counter(species2[gene_col])
keep_genes_s1 = []
keep_genes_s2 = []
for key in num_var_s1.keys():
    if num_var_s1[key] >= num_var_cut:
        keep_genes_s1.append(key)
        
for key in num_var_s2.keys():
    if num_var_s2[key] >= num_var_cut:
        keep_genes_s2.append(key)

if operation == "or":
    all_genes = list(np.unique(keep_genes_s1 + keep_genes_s2))
elif operation == "and":
    all_genes = list(np.unique(np.intersect1d(keep_genes_s1, keep_genes_s2)))
species1 = species1[species1[gene_col].isin(all_genes)]
species2 = species2[species2[gene_col].isin(all_genes)]

print(len(all_genes))

d_comp = {"A":"T", "C":"G", "G":"C", "T":"A"}
def revcomp(s):
    new_s = ""
    for i in s[::-1]:
        new_s = new_s + d_comp[i]
    return new_s
        
#Make trinuc equivalence classes
tris = list(set(species1[anc_tri_col]))
tris.sort()

d_equiv = {}
for i in tris:
    if i in d_equiv.values():
        d_equiv[i] = i
    else:
        d_equiv[i] = revcomp(i)

print(d_equiv)

def equiv(s):
    return d_equiv[s]
    
for i in d_equiv.keys():
    assert(d_equiv[i] == i or d_equiv[i] == revcomp(i))

species1[tri_col] = species1[tri_col].apply(equiv)
species1[anc_tri_col] = species1[anc_tri_col].apply(equiv)
species2[tri_col] = species2[tri_col].apply(equiv)

print(species1)

tris = np.unique(species1[anc_tri_col])

def gene_just_tri():
    ind = 1
    df_out = 0
    for i in list(np.unique(list(d_equiv.values()))):
        back = open(background_file)
        j = revcomp(i)
        assert(i == d_equiv[j])
        print(i, j)
        to_df = []
        if "mono" in control:
            for line in back:
                l = line.replace("\n", "").split("\t")
                if l[tri_col].upper()[1] == i or l[tri_col].upper()[1] == j:
                    to_df.append(l)
        else:
            for line in back:
                l = line.replace("\n", "").split("\t")
                if l[tri_col].upper() == i or l[tri_col].upper() == j:
                    to_df.append(l)
        #Count up the number of trinucs per gene
        #Subtract 1 to remove the added variant that ensured all genes would be present
        df_back = pd.DataFrame(to_df)
        df_back = df_back[df_back[gene_col].isin(all_genes)]
        counts_per_gene = Counter(list(df_back[gene_col]) + all_genes)
        counts_per_gene.subtract(Counter(all_genes))
        df2 = pd.DataFrame.from_dict(counts_per_gene, orient='index').reset_index().sort_values("index")
        df2.columns = ["Gene", "Number variants"]
        df2[i + "_" + j] = df2["Number variants"]/np.sum(df2["Number variants"])
        df2 = df2.sort_values("Gene").set_index("Gene")[[i + "_" + j]]
        if ind:
            ind = 0
            df_out = df2
        else:
            df_out = df_out.join(df2)
        back.close()
    if "mono" in control:
        df_out.to_csv(folder + "/" + "Species_MonoProbs_Back_All.txt", sep = "\t")
    else:
        df_out.to_csv(folder + "/" + "Species_TriProbs_Back_All.txt", sep = "\t")

if control == "just_tri" or control == "just_mono":
    out = []
    for i in tris:
        species1_tri = species1[species1[anc_tri_col].isin([i])].shape[0]
        species2_tri = species2[species2[tri_col].isin([i])].shape[0]
        out.append(species1_tri/(species1_tri + species2_tri))
    df = pd.DataFrame(out)
    df.index = tris
    df.columns = ["All_Scores"]
    if "mono" in control:
        df.to_csv(folder + "/Species_MonoProbs.txt", sep = "\t")
    else:
        df.to_csv(folder + "/Species_TriProbs.txt", sep = "\t")
 
    gene_just_tri()

    
else:
    #Get the bin boundaries using species1 + species2
    hist_bins = np.histogram(list(species1[phylop_col]) + list(species2[phylop_col]), bins = 100)[1]
    hist_bins[-1] = max_phylo
    print(list(hist_bins))
    
    #First pass to assign variant its PhyloP bin
    s1 = open(species1_file)
    out_s1 = open(species1_file.replace(".bed", "_Binned.bed"), 'w')
    
    for line in s1:
        l = line.replace("\n", "").split("\t")
        bin = np.digitize(np.float64(l[phylop_col]), hist_bins) - 1
        out_s1.write("\t".join(l + [str(bin)]) + "\n")
    s1.close()
    out_s1.close()
    
    s2 = open(species2_file)
    out_s2 = open(species2_file.replace(".bed", "_Binned.bed"), 'w')
    
    for line in s2:
        l = line.replace("\n", "").split("\t")
        bin = np.digitize(np.float64(l[phylop_col]), hist_bins) - 1
        out_s2.write("\t".join(l + [str(bin)]) + "\n")
    s2.close()
    out_s2.close()
    
    if control == "tri_both_score" or control == "mono_both_score":
        back = open(background_file)
        out_back = open(background_file.replace(".bed", "_Binned.bed"), 'w')
        for line in back:
            l = line.replace("\n", "").split("\t")
            bin = np.digitize(np.float64(l[phylop_col]), hist_bins) - 1
            out_back.write("\t".join(l + [str(bin)]) + "\n")
        back.close()
        out_back.close()
    
    out = []
    for i in tris:
        species1_tri = species1[species1[anc_tri_col].isin([i])]
        species2_tri = species2[species2[tri_col].isin([i])]
        species1_hist = np.histogram(species1_tri[phylop_col], bins = hist_bins)
        #Make a histogram with the same bin edges for species2
        species2_hist = np.histogram(species2_tri[phylop_col], bins = hist_bins)
        #Compute probability in each bin of the histogram
        #When there are zero total mutations, it will output zero
        #But that is fine as that will never be encountered later on
        species_prob_hist = species1_hist[0]/(np.maximum(species1_hist[0] + species2_hist[0], 1))
        out.append(species_prob_hist)
    df = pd.DataFrame(out)
    df.index = tris
    df.columns = hist_bins[0:len(hist_bins)-1]
    if "mono" in control:
        df.to_csv(folder + "/Species_MonoProbs.txt", sep = "\t")
    else:
        df.to_csv(folder + "/Species_TriProbs.txt", sep = "\t")
    
    if control == "tri_both_score" or control == "mono_both_score":
        for i in list(np.unique(list(d_equiv.values()))):
            
            to_df = []
            back = open(background_file)
            j = revcomp(i)
            assert(i == d_equiv[j])
            print(i, j)
            if "mono" in control:
                for line in back:
                    l = line.replace("\n", "").split("\t")
                    if l[tri_col].upper()[1] == i or l[tri_col].upper()[1] == j:
                        to_df.append(l)
            else:
                for line in back:
                    l = line.replace("\n", "").split("\t")
                    if l[tri_col].upper() == i or l[tri_col].upper() == j:
                        to_df.append(l)
            df_back = pd.DataFrame(to_df)
            df_back = df_back.sort_values(gene_col)
            df_back[phylop_col] = df_back[phylop_col].astype(np.float64)
            df_back = df_back[df_back[phylop_col] <= max_phylo]
            to_df = 0
           
            df2 = pd.DataFrame(df_back.groupby(gene_col).agg({phylop_col: lambda x: np.histogram(x, hist_bins)}))
            genes = list(df2.index)
            to_df2 = []
            [to_df2.append(x[0]) for x in list(df2[phylop_col])]
            df2 = pd.DataFrame(to_df2)
            df2.index = genes
            df2.columns = hist_bins[0:len(hist_bins)-1]
            
            #Adding this to restrict to only genes of interest and add back any missing genes
            df2 = df2.loc[np.intersect1d(df2.index, all_genes)]
            for gene in list(np.setdiff1d(all_genes, df2.index)):
                df2.loc[gene] = np.repeat(0, df2.shape[1])
            
            df2 = df2.sort_index()
            df2 = df2.divide(df2.sum(axis = 0), axis = 1)
            df2 = df2.replace(np.nan, np.float64(1/df2.shape[0]))
            if "mono" in control:
                df2.to_csv(folder + "/" + "Species_MonoProbs_Back_" + i + "_" + j + ".txt", sep = "\t")
            else:
                df2.to_csv(folder + "/" + "Species_TriProbs_Back_" + i + "_" + j + ".txt", sep = "\t")
            to_df2 = 0
            back.close()
    else:
        gene_just_tri()
        
                