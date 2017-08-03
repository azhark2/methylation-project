import pandas as pd
import pickle
import os
import pybedtools

# add 'gene_name' and 'gene description' column to each mutation file
#combine all mutation files into one file

df2 = pd.read_csv('gene_id_names.tsv', sep='\t')

#create dictionary that maps each Ensembl id to a tuple with (gene name, gene description)
idToGene = {}
for ids, name, description in zip(df2['Gene stable ID'], df2['Gene name'], df2['Gene description']):
    idToGene[ids] = (name, description)

files = [f for f in os.listdir('.') if os.path.isfile(f)]

for f in files:
    df = pd.read_csv(f, sep='\t', header=None)








