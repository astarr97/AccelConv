import pandas as pd
import numpy as np
from collections import Counter
from numpy.random import choice
from scipy.stats import mannwhitneyu as mwu
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
import sys
import os

species1_file = sys.argv[1]
species2_file = sys.argv[2]
start_iter = int(sys.argv[3])
end_iter = int(sys.argv[4])
write_perm = sys.argv[5]
name = sys.argv[6]

#Maximum allowable PhyloP score for non-background
#Typically 8.5 for our data
#Higher values indicate misalignment and/or other funny business
#Even though A/T nucleotides cannot generally have PhyloP scores this big, A/T nucleotides with PhyloP scores greater than 7.43 (the equivalent cutoff) do not necessarily indicate problems
#This is because they more likely indicate a nucleotide that was ancestrally A/T but became G/C in the common ancestor of the two species
max_phylo = float(sys.argv[7])

#What to control for during the simulation when assigning variants to species and gene
#Can be one of three options:
#(1) just_tri means only control for the trinucleotide context
#(2) tri_species_score means additionally controlling for the score (often phylop, but not necessarily) only when assigning to species
#(3) tri_both_score means additionally controlling for the score when assigning to both species and gene (not possible unless score is computed for all background variants)
#(4) just_mono means only control for the mononucleotide context (typically used for more closely related species)
#(5) mono_species_score is analogous to (2)
#(6) mono_both_score is analogous to (3)
#MUST BE THE SAME AS THE ARGUMENT TO make_prob_matrices.py
control = sys.argv[8]

#Folder containing the probability matrices to be used
prob_mat_fold = sys.argv[9]

#Minimum number of variants needed to be included
#Has format integer,and or integer,or; and implies both species must have at least integer variants; or implies just one species does
int_op = sys.argv[10]
num_var_cut = int(int_op.split(",")[0])
operation = int_op.split(",")[1]

#Whether we are doing this at the gene set level
#Input should be of the form file,maximum number of genes per category,minimum number of genes per category
try:
    gene_set_string = sys.argv[11]
    gene_set = gene_set_string.split(",")[0]
    min_genes = int(gene_set_string.split(",")[1])
    max_genes = int(gene_set_string.split(",")[2])
except:
    gene_set = 0

if write_perm == "0":
    write_perm = int(write_perm)
else:
    if write_perm not in os.listdir():
        os.mkdir(write_perm)

phylop_col = 1
gene_col = 2
tri_col = 3
anc_tri_col = 4

species1 = pd.read_csv(species1_file, sep = "\t", header = None).drop_duplicates()
species2 = pd.read_csv(species2_file, sep = "\t", header = None).drop_duplicates()

#Create a list that lets us remove anything containing an N (these overlapped with CDS for non-coding)
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
species1[anc_tri_col] = species1[anc_tri_col].str.upper()
species2[tri_col] = species2[tri_col].str.upper()

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
species1[phylop_col] = species1[phylop_col].astype(np.float64)
species2[phylop_col] = species2[phylop_col].astype(np.float64)

species1 = species1.drop(tri_col, axis = 1)
species1.columns = list(range(species1.shape[1]))

#Reset anc_tri_col as we have dropped tri_col
anc_tri_col = tri_col

d_comp = {"A":"T", "C":"G", "G":"C", "T":"A"}
def revcomp(s):
    new_s = ""
    for i in s[::-1]:
        new_s = new_s + d_comp[i]
    return new_s
        
#Make trinuc equivalence classes
tris = np.unique(species1[anc_tri_col])
tris.sort()

d_equiv = {}
for i in tris:
    if i in d_equiv.values():
        d_equiv[i] = i
    else:
        d_equiv[i] = revcomp(i)

def equiv(s):
    return d_equiv[s]
    
for i in d_equiv.keys():
    assert(d_equiv[i] == i or d_equiv[i] == revcomp(i))

species1[anc_tri_col] = species1[anc_tri_col].apply(equiv)
species2[tri_col] = species2[tri_col].apply(equiv)

#Usually superfluous, but we don't want negative PhyloP to count against genes
species1[phylop_col] = np.maximum(species1[phylop_col], 0)
species2[phylop_col] = np.maximum(species2[phylop_col], 0)

bin_index_col = 4

#Remove PhyloP scores that we know are wrong
species1 = species1[species1[phylop_col] <= max_phylo]
species2 = species2[species2[phylop_col] <= max_phylo]

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
print(species1)

#Also make sure we don't have any wrong binned variants
if control == "tri_both_score" or control == "tri_species_score" or control == "mono_both_score" or control == "mono_species_score":
    species1 = species1[species1[bin_index_col] < 100]
    species2 = species2[species2[bin_index_col] < 100]

species1["Total_Vars"] = np.repeat(1, species1.shape[0])
species2["Total_Vars"] = np.repeat(1, species2.shape[0])
species1_sum = species1.groupby([gene_col]).sum(numeric_only=1)
species2_sum = species2.groupby([gene_col]).sum(numeric_only=1)
if control == "just_tri" or control == "just_mono":
    species1_sum.columns = ["Species1 Sum PhyloP", "Num Species1 Var"]
    species2_sum.columns = ["Species2 Sum PhyloP", "Num Species2 Var"]
else:
    species1_sum.columns = ["Species1 Sum PhyloP", "Discard", "Num Species1 Var"]
    species2_sum.columns = ["Species2 Sum PhyloP", "Discard", "Num Species2 Var"]
    species1_sum = species1_sum.drop("Discard", axis = 1)
    species2_sum = species2_sum.drop("Discard", axis = 1)

new_sum = species1_sum.join(species2_sum, how = "outer").fillna(0)

def sum_gs(genes, df_metric):
    new_sum_genes = df_metric.loc[np.intersect1d(df_metric.index, genes.split(";"))]
    new_sum_genes = new_sum_genes.sum(axis = 0)
    return new_sum_genes

#If it is a gene set, do it at the gene set level
if gene_set:
    gs_df = pd.read_csv(gene_set, sep = "\t")
    keep_terms = []
    for index, row in gs_df.iterrows():
        genes_oi = np.intersect1d(row["Genes"].split(";"), new_sum.index)
        num_genes = len(genes_oi)
        if num_genes >= min_genes and num_genes <= max_genes:
            keep_terms.append([row["Term"], ";".join(list(genes_oi))])
    gs_df = pd.DataFrame(keep_terms)
    gs_df.columns = ["Term", "Genes"]

    new_sum_gs = gs_df["Genes"].apply(sum_gs, args=(new_sum, ))
    new_sum_gs.index = gs_df["Term"]
    new_sum_gs["PhyloP Difference"] = new_sum_gs["Species1 Sum PhyloP"] - new_sum_gs["Species2 Sum PhyloP"]
    new_sum_gs["PhyloP L2FC"] = np.log2((new_sum_gs["Species1 Sum PhyloP"] + 10)/(new_sum_gs["Species2 Sum PhyloP"] + 10))
    new_sum_gs.sort_values("PhyloP Difference")
    new_sum_gs.to_csv("AccelEvol_" + name + ".txt", sep = "\t")
    new_sum_gs["Total_Sims"] = np.repeat(0, new_sum_gs.shape[0])
    print(new_sum_gs)
else:
    new_sum["PhyloP Difference"] = new_sum["Species1 Sum PhyloP"] - new_sum["Species2 Sum PhyloP"]
    new_sum["PhyloP L2FC"] = np.log2((new_sum["Species1 Sum PhyloP"] + 10)/(new_sum["Species2 Sum PhyloP"] + 10))
    new_sum.sort_values("PhyloP Difference")
    new_sum.to_csv("AccelEvol_" + name + ".txt", sep = "\t")
    new_sum["Total_Sims"] = np.repeat(0, new_sum.shape[0])



#Read in the background probabilities
if "mono" in control:
    name_mat_add = "Mono"
else:
    name_mat_add = "Tri"
if control == "just_tri" or control == "tri_species_score" or control == "just_mono" or control == "mono_species_score":
    species_probs = pd.read_csv(prob_mat_fold + "/Species_" + name_mat_add + "Probs.txt", sep = "\t").set_index("Unnamed: 0")
    gene_probs = pd.read_csv(prob_mat_fold + "/Species_" + name_mat_add + "Probs_Back_All.txt", sep = "\t").set_index("Gene")
    if control == "tri_species_score":
        species_probs.columns = list(range(species_probs.shape[1]))
elif control == "tri_both_score" or control == "mono_both_score":
    species_probs = pd.read_csv(prob_mat_fold + "/Species_" + name_mat_add + "Probs.txt", sep = "\t").set_index("Unnamed: 0")
    d_gene_probs = {}
    for file in os.listdir(prob_mat_fold):
        if file != "Species_" + name_mat_add + "Probs.txt":
            v = pd.read_csv(prob_mat_fold + "/" + file, sep = "\t").set_index("Unnamed: 0")
            v.columns = list(range(v.shape[1]))
            d_gene_probs[file.split("_")[-2]] = v.copy()

#Commented out mann-whitney bit for run time while testing
"""to_join = []
for gene in list(new_sum.index):
    s1_gene = species1[species1[gene_col].isin([gene])]
    s2_gene = species2[species2[gene_col].isin([gene])]
    to_join.append([gene, np.median(s1_gene[phylop_col]), np.median(s2_gene[phylop_col]), mwu(s1_gene[phylop_col], s2_gene[phylop_col])[1]])
to_join = pd.DataFrame(to_join)
to_join.columns = ["Gene", "Median_PhyloP_Species1", "Median_PhyloP_Species2", "MWU p-value"]
new_sum = new_sum.join(to_join.set_index("Gene"))"""

#new_sum["Better_Sims_L2FC_wReal"] = np.repeat(0, new_sum.shape[0])
#new_sum["Better_Sims_Sum_wReal"] = np.repeat(0, new_sum.shape[0])
#new_sum["Better_Sims_L2FC_wSim"] = np.repeat(0, new_sum.shape[0])
#new_sum["Better_Sims_Sum_wSim"] = np.repeat(0, new_sum.shape[0])

#Important to center the lfc for things downstream to work, at least for now.
#new_sum["Centered PhyloP L2FC"] = new_sum["PhyloP L2FC"] - np.mean(new_sum["PhyloP L2FC"])
#new_sum["Centered PhyloP Difference"] = new_sum["PhyloP Difference"] - np.mean(new_sum["PhyloP Difference"])

if control == "just_tri" or control == "tri_species_score" or control == "just_mono" or control == "mono_species_score":
    assert(list(new_sum.index) == list(gene_probs.index))
else:
    for key in d_gene_probs.keys():
        assert(list(new_sum.index) == list(d_gene_probs[key].index))

#Prepare for permutations
input_sim = pd.concat([species1, species2]).sort_values(tri_col)
tris.sort()

np.random.seed(start_iter)

#Old functions for slower version commented out below (NOT USED)
def assign_to_species(x, prob_list):
    s1prob = prob_list[x]
    return choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob])
    
def assign_to_gene(x, gene_probs):
    prob_list = np.array(gene_probs.iloc[x])
    return choice(np.array(gene_probs.columns), p = prob_list)

print(input_sim)
tris = np.unique(input_sim[tri_col])
#For each permutation

for j in range(start_iter, end_iter):
    #print(j)
    #Go through the sorted list of trinucleotides and assign each mutation to a new gene and new species
    new_gene = []
    new_species = []
    for tri in tris:
        input_sim_tri = input_sim[input_sim[tri_col].isin([tri])].copy()
        #Slower version, no longer used!
        #if control == "just_tri" or control == "tri_species_score":
        #    gene_prob_dist = np.array(gene_probs[tri + "_" + revcomp(tri)])
        #    new_gene = new_gene + list(choice(np.array(gene_probs.index), p=gene_prob_dist, size = len(list(input_sim_tri.index))))
        #    if control == "just_tri":
        #        s1prob = species_probs.loc[tri]["All_Scores"]
        #        new_species = new_species + list(choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob], size = len(list(input_sim_tri.index))))
        #    elif control == "tri_species_score":
        #        prob_list = np.array(species_probs.loc[tri])
        #        new_species = new_species + list(input_sim_tri[bin_index_col].apply(assign_to_species, args = (prob_list, )))
        #elif control == "tri_both_score":
        #    prob_list = np.array(species_probs.loc[tri])
        #    new_species = new_species + list(input_sim_tri[bin_index_col].apply(assign_to_species, args = (prob_list, )))
        #    gene_probs = d_gene_probs[tri].T
        #    new_gene = new_gene + list(input_sim_tri[bin_index_col].apply(assign_to_gene, args = (gene_probs, )))
        #    print(tri)
        #input_sim["New_Gene"] = new_gene
        #input_sim["New_Species"] = new_species
        #print(len(new_gene), len(new_species), len(tris))
        #print(input_sim)
        
        #Newer version, much faster!
        if control == "just_tri" or control == "just_mono":
            gene_prob_dist = np.array(gene_probs[tri + "_" + revcomp(tri)])
            input_sim_tri["New_Gene"] = list(choice(np.array(gene_probs.index), p=gene_prob_dist, size = len(list(input_sim_tri.index))))
            s1prob = species_probs.loc[tri]["All_Scores"]
            input_sim_tri["New_Species"] = list(choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob], size = len(list(input_sim_tri.index))))
        elif control == "tri_species_score" or control == "mono_species_score":
            #Sort things by bin_index
            input_sim_tri = input_sim_tri.sort_values(bin_index_col)
            bin_indices = np.unique(input_sim_tri[bin_index_col])
            
            #Identical to just_tri for New_Gene
            gene_prob_dist = np.array(gene_probs[tri + "_" + revcomp(tri)])
            input_sim_tri["New_Gene"] = list(choice(np.array(gene_probs.index), p=gene_prob_dist, size = len(list(input_sim_tri.index))))
            
            #Iterate through the bin_indices
            new_species = []
            prob_list = np.array(species_probs.loc[tri])
            for bin_index in bin_indices:
                #Get the size of the sample (number of species differences with the bin_index and trinucleotide)
                size_to_sample = len(list(input_sim_tri[input_sim_tri[bin_index_col].isin([bin_index])].index))
                s1prob = prob_list[bin_index]
                new_species = new_species + list(choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob], size = size_to_sample))
            input_sim_tri["New_Species"] = new_species
        elif control == "tri_both_score" or control == "mono_both_score":
            #Sort things by bin_index
            input_sim_tri = input_sim_tri.sort_values(bin_index_col)
            bin_indices = np.unique(input_sim_tri[bin_index_col])

            #Iterate through the bin_indices and assign genes
            new_species = []
            new_gene = []
            prob_list = np.array(species_probs.loc[tri])
            gene_probs = d_gene_probs[tri].T

            for bin_index in bin_indices:
                size_to_sample = len(list(input_sim_tri[input_sim_tri[bin_index_col].isin([bin_index])].index))
                gene_prob_list = np.array(gene_probs.iloc[bin_index])
                
                s1prob = prob_list[bin_index]
                new_gene = new_gene + list(choice(np.array(gene_probs.columns), p = gene_prob_list, size = size_to_sample))
                new_species = new_species + list(choice(["S1Der", "S2Der"], p = [s1prob, 1-s1prob], size = size_to_sample))
            input_sim_tri["New_Gene"] = new_gene
            input_sim_tri["New_Species"] = new_species
            #print(tri)
        if tri == tris[0]:
            input_sim_new = input_sim_tri
        else:
            input_sim_new = pd.concat([input_sim_new, input_sim_tri])
    
    #print(input_sim_new)
    #Checking it worked for just_tri
    #s1tot = input_sim_new[(input_sim_new[tri_col].isin(["TTT"])) & input_sim_new["New_Species"].isin(["S1Der"])].shape[0]
    #s2tot = input_sim_new[(input_sim_new[tri_col].isin(["TTT"])) & input_sim_new["New_Species"].isin(["S2Der"])].shape[0]
    #print(s1tot/(s2tot + s1tot), s1tot, s2tot)
    
    input_sim_species1 = input_sim_new[input_sim_new["New_Species"].isin(["S1Der"])]
    input_sim_species2 = input_sim_new[input_sim_new["New_Species"].isin(["S2Der"])]
    #input_sim_species1.to_csv("Test_Sim_S1_" + name + ".bed", sep = "\t", header = False, index = False)
    #input_sim_species2.to_csv("Test_Sim_S2_" + name + ".bed", sep = "\t", header = False, index = False)

    #Do the same computation as done for the original
    species1_sum = input_sim_species1.groupby(["New_Gene"]).sum(numeric_only=1)
    species2_sum = input_sim_species2.groupby(["New_Gene"]).sum(numeric_only=1)
    
    if control == "just_tri" or control == "just_mono":
        species1_sum.columns = ["Species1 Sum PhyloP Sim", "Num Species1 Var Sim"]
        species2_sum.columns = ["Species2 Sum PhyloP Sim", "Num Species2 Var Sim"]
    else:
        species1_sum.columns = ["Species1 Sum PhyloP Sim", "Discard", "Num Species1 Var Sim"]
        species2_sum.columns = ["Species2 Sum PhyloP Sim", "Discard", "Num Species2 Var Sim"]
        species1_sum = species1_sum.drop("Discard", axis = 1)
        species2_sum = species2_sum.drop("Discard", axis = 1)
    
    new_sim = species1_sum.join(species2_sum, how = "outer").fillna(0)
    if gene_set:
        new_sim_gs = gs_df["Genes"].apply(sum_gs, args=(new_sim, ))
        new_sim_gs.index = gs_df["Term"]

        new_sum_gs = new_sum_gs.join(new_sim_gs)
        new_sum_gs = new_sum_gs.fillna(0)
        new_sum_gs["PhyloP Difference Sim"] = new_sum_gs["Species1 Sum PhyloP Sim"] - new_sum_gs["Species2 Sum PhyloP Sim"]
        new_sum_gs["PhyloP L2FC Sim"] = np.log2((new_sum_gs["Species1 Sum PhyloP Sim"] + 10)/(new_sum_gs["Species2 Sum PhyloP Sim"] + 10))
    
        #If desired, write out the new_sim because we are going to be using the normal approximation for this one
        if write_perm:
            new_sum_gs.to_csv(write_perm + "/Simulation" + str(j) + "_" + name + ".txt", sep = "\t")

        new_sum_gs["Total_Sims"] = new_sum_gs["Total_Sims"] + 1
        new_sum_gs = new_sum_gs.drop(["Species1 Sum PhyloP Sim", "Num Species1 Var Sim", "Species2 Sum PhyloP Sim", "Num Species2 Var Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim"], axis = 1)

    else:
        #Join together, anything that didn't have a value in the simulation was assigned no mutations
        new_sum = new_sum.join(new_sim)
        new_sum = new_sum.fillna(0)
        new_sum["PhyloP Difference Sim"] = new_sum["Species1 Sum PhyloP Sim"] - new_sum["Species2 Sum PhyloP Sim"]
        new_sum["PhyloP L2FC Sim"] = np.log2((new_sum["Species1 Sum PhyloP Sim"] + 10)/(new_sum["Species2 Sum PhyloP Sim"] + 10))
    
        #If desired, write out the new_sim because we are going to be using the normal approximation for this one
        if write_perm:
            new_sum.to_csv(write_perm + "/Simulation" + str(j) + "_" + name + ".txt", sep = "\t")

        new_sum["Total_Sims"] = new_sum["Total_Sims"] + 1
        new_sum = new_sum.drop(["Species1 Sum PhyloP Sim", "Num Species1 Var Sim", "Species2 Sum PhyloP Sim", "Num Species2 Var Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim"], axis = 1)

new_sum.to_csv("AccelEvol_" + name + "_" + str(start_iter) + "-" + str(end_iter) + ".txt", sep = "\t")


### OLD STUFF COMMENTED OUT BEFORE USING NORMAL APPROX FOR EVERYTHING ###
### WOULD NEED TO MOVE UP ###
#Add 1 for each time the permutation had a more extreme score than the real value

#new_sum["Centered wReal PhyloP L2FC Sim"] = new_sum["PhyloP L2FC Sim"] - np.mean(new_sum["PhyloP L2FC"])
#new_sum["Centered wReal PhyloP Difference Sim"] = new_sum["PhyloP Difference Sim"] - np.mean(new_sum["PhyloP Difference"])
#new_sum["Centered wSim PhyloP L2FC Sim"] = new_sum["PhyloP L2FC Sim"] - np.mean(new_sum["PhyloP L2FC Sim"])
#new_sum["Centered wSim PhyloP Difference Sim"] = new_sum["PhyloP Difference Sim"] - np.mean(new_sum["PhyloP Difference Sim"])

#new_sum["Better_Sims_L2FC_wReal"] = new_sum["Better_Sims_L2FC_wReal"] + pd.DataFrame(((new_sum["Centered PhyloP L2FC"] <= new_sum["Centered wReal PhyloP L2FC Sim"]) & (new_sum["Centered PhyloP L2FC"] >= 0)) | ((new_sum["Centered PhyloP L2FC"] >= new_sum["Centered wReal PhyloP L2FC Sim"]) & (new_sum["Centered PhyloP L2FC"] < 0))).replace({True:1, False:0})[0]
#new_sum["Better_Sims_Sum_wReal"] = new_sum["Better_Sims_Sum_wReal"] + pd.DataFrame(((new_sum["Centered PhyloP Difference"] <= new_sum["Centered wReal PhyloP Difference Sim"]) & (new_sum["Centered PhyloP Difference"] >= 0)) | ((new_sum["Centered PhyloP Difference"] >= new_sum["Centered wReal PhyloP Difference Sim"]) & (new_sum["Centered PhyloP Difference"] < 0))).replace({True:1, False:0})[0]
#new_sum["Better_Sims_L2FC_wSim"] = new_sum["Better_Sims_L2FC_wSim"] + pd.DataFrame(((new_sum["Centered PhyloP L2FC"] <= new_sum["Centered wSim PhyloP L2FC Sim"]) & (new_sum["Centered PhyloP L2FC"] >= 0)) | ((new_sum["Centered PhyloP L2FC"] >= new_sum["Centered wSim PhyloP L2FC Sim"]) & (new_sum["Centered PhyloP L2FC"] < 0))).replace({True:1, False:0})[0]
#new_sum["Better_Sims_Sum_wSim"] = new_sum["Better_Sims_Sum_wSim"] + pd.DataFrame(((new_sum["Centered PhyloP Difference"] <= new_sum["Centered wSim PhyloP Difference Sim"]) & (new_sum["Centered PhyloP Difference"] >= 0)) | ((new_sum["Centered PhyloP Difference"] >= new_sum["Centered wSim PhyloP Difference Sim"]) & (new_sum["Centered PhyloP Difference"] < 0))).replace({True:1, False:0})[0]

#new_sum["Better_Sims_L2FC_wReal"] = new_sum["Better_Sims_L2FC_wReal"] + pd.DataFrame(np.abs(new_sum["Centered PhyloP L2FC"]) <= np.abs(new_sum["Centered wReal PhyloP L2FC Sim"])).replace({True:1, False:0})[0]
#new_sum["Better_Sims_Sum_wReal"] = new_sum["Better_Sims_Sum_wReal"] + pd.DataFrame(np.abs(new_sum["Centered PhyloP Difference"]) <= np.abs(new_sum["Centered wReal PhyloP Difference Sim"])).replace({True:1, False:0})[0]
#new_sum["Better_Sims_L2FC_wSim"] = new_sum["Better_Sims_L2FC_wSim"] + pd.DataFrame(np.abs(new_sum["Centered PhyloP L2FC"]) <= np.abs(new_sum["Centered wSim PhyloP L2FC Sim"])).replace({True:1, False:0})[0]
#new_sum["Better_Sims_Sum_wSim"] = new_sum["Better_Sims_Sum_wSim"] + pd.DataFrame(np.abs(new_sum["Centered PhyloP Difference"]) <= np.abs(new_sum["Centered wSim PhyloP Difference Sim"])).replace({True:1, False:0})[0]

#new_sum = new_sum.drop(["Species1 Sum PhyloP Sim", "Num Species1 Var Sim", "Species2 Sum PhyloP Sim", "Num Species2 Var Sim", "PhyloP Difference Sim", "PhyloP L2FC Sim", "Centered wReal PhyloP Difference Sim", "Centered wReal PhyloP L2FC Sim", "Centered wSim PhyloP Difference Sim", "Centered wSim PhyloP L2FC Sim"], axis = 1)
