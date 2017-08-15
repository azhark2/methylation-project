
import pandas as pd
import pickle
import os
import csv

# add 'gene_name' and 'gene description' column to each mutation file
#combine all mutation files into one file


#create dictionary that maps each Ensembl id to a tuple with (gene name, gene description)

idToGene = pickle.load(open('idToGene.pickle', 'rb'))
os.chdir('/data/khandekara2/bed_CpGs/mutations')
files = [f for f in os.listdir('/data/khandekara2/bed_CpGs/mutations') if os.path.isfile(f)]


with open('all_mutations.bed', 'w') as csvout:
    writer = csv.writer(csvout, delimiter='\t')
    for file in files:
           with open(file, 'r') as f:
               reader = csv.reader(f, delimiter='\t')
               for row in reader:
                  if row[7] in idToGene.keys():
                     gene_name = idToGene[row[7]][0]
                     gene_function = idToGene[row[7]][1]
                     row.append(gene_name)
                     row.append(gene_function)
                  else:
                     row.append('MISSING')
                     row.append('MISSING')

                  writer.writerow(row)








