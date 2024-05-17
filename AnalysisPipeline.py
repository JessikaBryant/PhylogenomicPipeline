###Pipeline to assemble genomes, infer gene trees and genome trees to investigate phylogenetic incongruence.
#load modules
import os
import sys
import glob
import csv
import matplotlib.pyplot as plt 
import subprocess
from Bio import SeqIO
from Bio import Phylo
from Bio.Seq import Seq

###Input Gene family sequences
indir ="/shared/forsythe/BB485/Week06/Brass_CDS_seqs/" #input directory where fasta files live
outdir="/scratch/bryantj2/bb485/week06/PhylogenomicPipeline/OutputFolder/" #output directory where future files will go 
filelist=glob.glob(indir + "*.fasta")#list of files
#print(len(filelist))



#####Perform Multiple Sequence Alignment (Mafft)#####
#for loop to create an output for each file
for file in filelist:
    #print(file)
    new_file_path = file.replace(indir, outdir) #make  filepath
    #print(new_file_path) #print file path
    #Create system call for mafft
    aln_cmd = "mafft --auto --quiet "+file +"> "  +new_file_path
    #Check the command 
    #print(aln_cmd)
    #Run the command 
    os.system(aln_cmd) #Uncomment this once you’ve double-checked that it’s looking good.



#####Infer Phylogeny of Gene Tree (IQ Tree)#####

#list Aligned Fasta files
aln_files=glob.glob(outdir + "*.fasta")
#print(aln_files)

for aln in aln_files:
    #Create the command. -nt 2 means two threads. If running this from within a job submission, you could use more threads to make it go faster.
    tree_command = f"iqtree -s {aln} -m TEST -nt 2"
    #Check the command 
    #print(tree_command)
    #Run the command using a 'system call'
    os.system(tree_command) #uncomment once you've check the command



#####Output a visual Gene Tree (Biopython)######
#list trees from output of above
TreeFiles=glob.glob(outdir+"*.treefile")
#print(TreeFiles)

#create place for topologies to go
topo12=0
topo13=0
topo23=0
Unknown=0   

#For loop to read in tree files and parse through them to get topologies of trees
for tree in TreeFiles:
    #Read in the tree and store as phylo object
    temp_tree = Phylo.read(tree, "newick")
    print(temp_tree)
    #Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip
            #Stope the loop once we found the correct tip
            break
        
    #Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
        
    #Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
        
    #Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp=t 
        elif "Cr_" in t.name:
            Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t
            
    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp] 

    #Use series of if/else statemetns to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
        topo12+=1
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
        topo13+=1
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
        topo23+=1
    else:
        topo_str = "Unknown"
        Unknown+=1

    fields=["12top", "13top", "23top","Unknown"]
    rows=str([topo12,topo13,topo23,Unknown])
    #write topo_str to an outfile     
    with open(outdir+"TreeSearch.csv", "w") as f:
        write = csv.writer(f)
        
        write.writerow(fields)
        write.writerows(rows)

    #print(topo13)
    #print(topo12)
    #print(topo23)
    #print(Unknown)

    #print(topo_str)


"""
####Create Pie Chart from topos####
#assign labels for pie chart
labels= "12top", "13top", "23top", "Unknown"
#print(labels)
#count and assign the percentages of each topo
sizes= [len("12top"), len("13top"), len("23top"), len("Unknown")]
#print(sizes)
#start plot
fig, ax = plt.subplots()
#create plot
ax.pie(sizes, labels=labels, autopct='%1.1f%%')
#save plot
plt.savefig("treechart.pdf", format='pdf')
"""