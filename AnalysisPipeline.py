###Pipeline to assemble genomes, infer gene trees and genome trees to investigate phylogenetic incongruence.
#load modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio import Phylo
from Bio.seq import Seq

###Input Gene family sequences
indir =""
outdir="/scratch/bryantj2/bb485/week06/PhylogenomicPipeline/OutputFolder"

###Perform Multiple Sequence Alignment (Mafft)

###Infer Phylogeny of Gene Tree (IQ Tree)

###Output a visual Gene Tree (Biopython)
