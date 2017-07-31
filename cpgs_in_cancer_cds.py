#counts the total # of CpG's in the cds of every Cosmic cancer gene

import pybedtools
import os
import csv

with open('cancer_genes_hg19.bed', 'r') as f, open('results.csv', 'w') as result:
    reader = csv.reader(f, delimiter='\t')
    writer = csv.writer(result, delimiter='\t')
    for row in reader:
        with open('temp.bed', 'w') as csvout:
            writer2 = csv.writer(csvout, delimiter='\t')
            writer2.writerow(row)
        gene_name = row[3]
        a = pybedtools.BedTool('all_cds_cpgs_shifted_2.bed', )
        a.intersect('temp.bed').saveas(gene_name + 'cds_cpgs.bed')
        with open(gene_name + 'cds_cpgs.bed', 'r') as file:
            writer.writerow([gene_name, len(file.readlines())])
