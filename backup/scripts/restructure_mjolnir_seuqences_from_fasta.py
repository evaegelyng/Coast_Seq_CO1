from pysam import FastaFile
import pandas as pd

X = pd.read_csv("./classified_mjolnir_csv.csv", sep="\t") #read table
sequences_object = FastaFile("COSQ_seeds_abundant.fasta") #read fasta file
to_add = [sequences_object.fetch(i) for i in X['id']] #read each sequence in fasta file from table IDs
X['sequence'] = to_add
X.to_csv("./COSQ_classified.tsv", sep="\t", index=None, na_rep='NA')