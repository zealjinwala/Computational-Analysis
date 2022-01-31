# By Zeal Jinwala
# Date: January 30, 2022
# Data Source: The Illumina and nanopore sequence datasets of the nose swab samples, generated and analyzed
#               in the current study, are available in the European Nucleotide Archive (ENA) under accession number
#               PRJEB28612. https://www.ebi.ac.uk/ena/browser/view/PRJEB28612?show=reads 

import pandas as pd
import matplotlib as plt


# input Nanopore Seq Data
# read input
def fasta_readfirst(file):
    from Bio import SeqIO
    df = pd.DataFrame(columns=['IDs', 'Seqs'])
    print("here")
    for record in SeqIO.parse(file, "fastq"):
        df = df.append({'IDs': record.id, 'Seqs': record.seq}, ignore_index=True)
    return df

fileIllumina = "/Users/zsj24/GitHub/Computational-Analysis/sampleIllumina.fastq"
fileNanopore = "/Users/zsj24/GitHub/Computational-Analysis/sampleNanopore.fastq"
dfIllumina = fasta_readfirst(fileIllumina)
dfNanopore = fasta_readfirst(fileNanopore)

# read length distribution
readLengthsIL = dfIllumina['Seqs'].str.len()
readLengthsNP = dfNanopore['Seqs'].str.len()
plt.pyplot.hist(readLengthsIL)
# quality score distribution 
# GC content Histogram

# input Illumina Seq Data
# read input
# read length distribution
# quality score distribution 
# GC content Histogram

# differences between Illumina and Nanopore 


# filter sequencing reads