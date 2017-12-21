#plots the methylation profiles at each site for all genes in the cancer gene census. The mean methylation ratio and standard deviation are shown for each site(aggregated across all samples for that site)


import matplotlib.pyplot as plt
import pybedtools
import pandas as pd
import csv
import fileinput

def plot_profile(file, cancer_type, gene):
    df = pd.read_csv(file, sep='\t')
    df = pd.read_csv(file, sep='\t')
    grouped = df.groupby('location')
    groupby = grouped.agg({'methylation_ratio': [np.mean, np.std], 'prev_base': 'max'})
    groupby.reset_index(inplace=True)
    groupby.columns = groupby.columns.droplevel(0)
    # groupby.to_csv(file[:-4] + '.tsv', sep='\t')
    y = groupby['mean']
    y_err = groupby['std']
    x = [i for i in range(1, groupby.shape[0] + 1)]
    x_labels = [i + 'CG' for i in groupby['max']]
    fig = plt.figure(figsize=(20, 2))
    ax = fig.add_subplot(111)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation='vertical', fontsize=12)
    ax.errorbar(x, y, yerr=y_err, fmt='o')
    ax.set_xlabel('Coordinate', fontsize=14)
    ax.set_ylabel('Average Methylation Ratio', fontsize=14)
    ax.scatter(x, y2, color='red')
    ax.set_title('Methylation Profile of TP53 Coding Region in MALY', fontsize=16)





a = pybedtools.BedTool('all_cds_cpgs_shifted.bed')
pbca = pybedtools.BedTool('PBCA_prevBase.bed')
maly = pybedtools.BedTool('MALY_prevBase.bed')

genesToStdDev = {}


PBCA_files = []
MALY_files = []
with open('cancer_genes.bed', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        with open('temp_file.bed', 'w') as csvout:
            writer = csv.writer(csvout, delimiter='\t')
            writer.writerow([row[0], row[1], row[2]])

            a.intersect('temp_file.bed').saveas(row[3] + '_cds_CpGs.bed')
            pbca.intersect(row[3] + '_cds_CpGs.bed').saveas('PBCA_' + row[3] + '_cds_CpGs.bed')
            maly.intersect(row[3] + '_cds_CpGs.bed').saveas('MALY_' + row[3] + '_cds_CpGs.bed')

            PBCA_files.append('PBCA_' + row[3] + '_cds_CpGs.bed')
            MALY_files.append('MALY_' + row[3] + '_cds_CpGs.bed')

#add headers
headers = 'chromosome start stop id methylation_ratio prev_base'.split()
for file in PBCA_files:
    for line in fileinput.input([file], inplace=True):
        print ('\t'.join(headers))
        print (line)


    plot_profile(file, 'PBCA_', file.split('_')[1])

for file in MALY_files:
    for line in fileinput.input([file], inplace=True):
        print('\t'.join(headers))
        print(line)
    plot_profile(file, "MALY_")











#list of cosmic cancer genes

