cd###Pipeline to assemble genomes, infer gene trees and genome trees to investigate phylogenetic incongruence.
#load modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio import Phylo
from Bio.Seq import Seq

###Input Gene family sequences
indir ="/shared/forsythe/BB485/Week06/Brass_CDS_seqs/" #input directory where fasta files live
outdir="/scratch/bryantj2/bb485/week06/PhylogenomicPipeline/OutputFolder" #output directory where future files will go
file ="" #empty file 

filelist=glob.glob(indir + "*.fasta")#list of files
#print(len(filelist))

###Perform Multiple Sequence Alignment (Mafft)
#for loop to create an output for each file
for filename in filelist:
    file ="" #empty file 
    new_file_path = file.replace(indir, outdir) #make  filepath
    print(new_file_path) #print file path
    #Create system call for mafft
    aln_cmd = "mafft –auto –quiet "+file+" > " +new_file_path
    #Check the command 
    print(aln_cmd)
    break
sys.exit()




#Run the command 
#os.system(aln_cmd) #Uncomment this once you’ve double-checked that it’s looking good.


###Infer Phylogeny of Gene Tree (IQ Tree)
#Create the command. -nt 2 means two threads. If running this from within a job submission, you could use more threads to make it go faster.
tree_command = f"iqtree -s {aln} -m TEST -nt 2"

#Check the command 
print(tree_command)


#Run the command using a 'system call'
# os.system(tree_command) #uncomment once you've check the command

###Output a visual Gene Tree (Biopython)
    
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
elif bool(temp_tree.is_monophyletic(P1_and_P3)):
    topo_str = "13top"
elif bool(temp_tree.is_monophyletic(P2_and_P3)):
    topo_str = "23top"
else:
    topo_str = "Unknown"

print(topo_str)