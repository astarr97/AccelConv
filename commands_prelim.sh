#Commands to generate the input bed files
#Finds SNPs, computes PhyloP scores, and intersects

#focal_species is the species that has the trait of interest
focal_species="Equus_caballus"
#comp_species is the closest living relative of the focal_species that does not have the trait of interest
comp_species="Equus_asinus"
#out_species is the out_group used to determine if a trait is focal_species-derived or comp_species-derived
out_species="Tapirus_terrestris"

#These are additional species that we might want to also find SNPs for
#Must be comma-separated list of species
#species_add_foc is other species that are more closely related to the focal_species than any other species and share the trait of interest
species_add_foc="Equus_przewalskii"
#species_add_rel is other species that are more closely related to the comp_species than the outgroup or the focal_species
species_add_rel="NA"
#species_add_out are other viable outgroups
species_add_out="Diceros_bicornis"

#For example, if you were interested in aquatic living, the focal_species might be Orcinus_orca, the comp_species might be Bos_taurus, and the out_species might be Sus_scrofa
#In this case, all cetaceans are in the same clade as Orcinus_orca and live under water, so could be species_add_foc
#All ruminants are in he same clade as Bos_taurus, so could be species_add_rel
#The species_add_out could be another Sus species or any species further from cetaceans/ruminants

#Create a fasta file for the focal species
/home/groups/hbfraser/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta ../hg38.447way.hal $focal_species > $focal_species.fasta
samtools faidx $focal_species.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $focal_species.fasta.fai | sort -n -r -k 3,3 > $focal_species.bed

#Command to generate the scripts to actually do the computations
cp /home/groups/hbfraser/astarr_scripts/AccelConv/make_scripts_conv.py ./
python make_scripts_conv.py $focal_species.bed $focal_species $comp_species $out_species $species_add_foc $species_add_rel $species_add_out

#Submit jobs to run the scripts
cp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh /scratch/users/astarr97/PhyloP/Feasibility/$focal_species
cd $focal_species
./driver.sh
cd ..
