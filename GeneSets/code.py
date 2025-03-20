import gseapy as gs
import numpy as np
import gseapy as gs

#This code was used to generate these files

HPO = pd.read_csv("../PosSelect_Testing/HPO_World_14_Typ.csv").set_index("HPO_IDs")
HPO_Conv = pd.read_csv("../PosSelect_Testing/HPO_Conversion.csv")
HPO_Conv = HPO_Conv.set_index("HPO_IDs")
HPO = HPO_Conv.join(HPO).set_index("HPO_Phenos")
out = []
for index, row in HPO.iterrows():
    out.append([index, ";".join(list(HPO.columns[np.where(row)]))])
    
df = pd.DataFrame(out)
df.columns = ["Term", "Genes"]
df.to_csv("HPO_AccelEvol_Input.txt", sep = "\t", index = False)

process = gs.get_library(name='GO_Biological_Process_2023', organism='Human')
out = open("GOBP_AccelEvol_Input.txt", 'w')
out.write("Term\tGenes\n")
for key in process.keys():
    out.write("\t".join([key, ";".join(process[key])]) + "\n")
out.close()

process = gs.get_library(name='GO_Molecular_Function_2023', organism='Human')
out = open("GOMF_AccelEvol_Input.txt", 'w')
out.write("Term\tGenes\n")
for key in process.keys():
    out.write("\t".join([key, ";".join(process[key])]) + "\n")
out.close()

process = gs.get_library(name='GO_Cellular_Component_2023', organism='Human')
out = open("GOCC_AccelEvol_Input.txt", 'w')
out.write("Term\tGenes\n")
for key in process.keys():
    out.write("\t".join([key, ";".join(process[key])]) + "\n")
out.close()

process = gs.get_library(name='KEGG_2021_Human', organism='Human')
out = open("KEGG_AccelEvol_Input.txt", 'w')
out.write("Term\tGenes\n")
for key in process.keys():
    out.write("\t".join([key, ";".join(process[key])]) + "\n")
out.close()

process = gs.get_library(name='CORUM', organism='Human')
out = open("CORUM_AccelEvol_Input.txt", 'w')
out.write("Term\tGenes\n")
for key in process.keys():
    out.write("\t".join([key, ";".join(process[key])]) + "\n")
out.close()

organizer = pd.read_csv("../PosSelect_Testing/ORGANizer_World_TYP_v14.csv").set_index("Gene").T
out = []
for index, row in organizer.iterrows():
    out.append([index, ";".join(list(HPO.columns[np.where(row)]))])
df = pd.DataFrame(out)
df.columns = ["Term", "Genes"]
df.to_csv("GeneOrganizer_AccelEvol_Input.txt", sep = "\t", index = False)

v = pd.read_csv("gnomad.v4.1.constraint_metrics.tsv", sep = "\t")
v = v[v["lof.pLI"] > 0.9]
keep_znf = []
for i in np.unique(v["gene"].dropna()):
    if "ZNF" in i:
        keep_znf.append(i)
pd.DataFrame(keep_znf).to_csv("Keep_ZNF_Genes.txt", sep = "\t", header = None, index = False)