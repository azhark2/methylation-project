import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.stats import binom_test
import os
import pybedtools
#first apply binomial test to each gene and write out results to file
path = '/data/khandekara2/normal_WGBS/processed_data/methylation_mutation/'

# for file in os.listdir(path):
#     if file.endswith('.cds.CpGs'):
#         tissue = file.split('.')[1]
#         a = pybedtools.BedTool(path + file)
#         df = pd.read_csv(path + file, sep='\t') # this step will take a long execution time
#         p = df.iloc[:, 4].mean() # p is the average methylation across genome
#         with open('/data/khandekara2/bed_CpGs/mart_export.bed', 'r') as csvin:
#             reader = csv.reader(csvin, delimiter='\t')
#             no_Cpgs = 0
#             for row in reader: #for every protein-coding gene
#                 gene_name = row[5]
#                 with open('current_gene.bed', 'w') as csvout:
#                     writer = csv.writer(csvout, delimiter='\t')
#                     writer.writerow(row)
#                 b = pybedtools.BedTool('current_gene.bed')
#                 a.intersect(b, wa=True, wb=True).saveas(tissue + '_' + gene_name + '_cds_CpGs.bed') # all CpG's in coding region of given gene


#create dictionary that maps protein coding gene to it's length(defined as translation start to stop codon)
geneToLength = {}
genes = pd.read_csv('/data/khandekara2/bed_CpGs/mart_export.bed', sep='\t', header=None)
for gene, length in zip(genes.iloc[:,3], genes.iloc[:,5]):
    geneToLength[gene] = length

cdsAvg = pickle.load(open('cdsAvg.pickle', 'rb')) #dictionary that maps benign tissue type to average methylation level in cds

with open('P_cg_results_highly_methylated.tsv', 'w') as output:
    writer1 = csv.writer(output, delimiter='\t')
    for file in os.listdir('.'):
        if file.endswith('.cds.CpGs'):
            tissue = file.split('_')[0]
            p = cdsAvg[tissue]
            gene_name = file.split('_')[1]
            if os.stat(tissue + '_' + gene_name + '_cds_CpGs.bed').st_size != 0:
                print(gene_name + 'had no CpGs')
                no_Cpgs += 1
                writer1.writerow([tissue, gene_name, 0, 0, 0, None, p, geneToLength[gene]])
            else:
                df = pd.read_csv(file, sep='\t', header=None)
                length = df.iloc[2, 9]
                methylated = df[df.iloc[:, 4] >= 0.8]
                unmethylated = df[df.iloc[:, 4] <= 0.2]
                fuzzy = df[df.iloc[:, 4] > 0.2 & df.iloc[:, 4] < 0.8]
                nCpG = unmethylated.shape[0] + methylated.shape[0]
                mCpG = methylated.shape[0]
                iCpG = fuzzy.shape[0]
                P_CpG = binom_test(mCpG, nCpG, p, alternative='greater')
                writer1.writerow([tissue, gene_name, nCpG, mCpG, iCpG, P_CpG, p, length]) #tissue, gene name, nCpG, mCpG, probability of observing this deviation by chance, average genome methylation


# #now plot results
# df = pd.read_csv('P_cg_results_highly_methylated.tsv', sep='\t', header=None)
# fig = plt.figure(figsize=(17,15), dpi=300)
# fig.suptitle('p is average methylation of cds, alternative = greater')
# for i, tissue in enumerate(df.iloc[:, 0].unique()): #plot for each tissue
#     sub_df = df[df.iloc[:, 0] == tissue]
#     ax = fig.add_subplot(4, 3, i+1)
#     p_values = np.array(df.iloc[:, 5])
#     n,bins,edges = ax.hist(p_values, bins=20, weights=np.zeros_like(p) + 1. / p_values.size)
#     ax.set_xticks([0.05, 0.95])
#     ax.set_ylabel('Frequency', fontsize=14)
#     ax.set_xlabel(r'$P_{CG}$', fontsize=18)
#
# fig.savefig('P_cg_results_highly_methylated.png')
