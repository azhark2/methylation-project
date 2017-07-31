#plots the methylation profiles at each site for all genes in the cancer gene census. The mean methylation ratio and standard deviation are shown for each site(aggregated across all samples for that site)


import matplotlib.pyplot as plt
import pybedtools
import pandas as pd
import csv
import fileinput

def plot_profile(file, cancer_type, gene):
    df = pd.read_csv(file, sep='\t')
    locations = []
    for chrom, start, stop in zip(list(df['chromosome']), list(df['start']), list(df['stop'])):
        locations.append((chrom, start, stop))
    df['location'] = locations

    grouped = df['methylation_ratio'].groupby(df['location'])
    groupby = grouped.agg([np.mean, np.std])
    groupby.reset_index(inplace=True)
    groupby.to_csv(file[:-4] + '.tsv', sep='\t')
    # plot
    y = groupby['mean']
    y_err = groupby['std']
    avg = sum(y_err) / float(len(y_err))
    genesToStdDev[gene] = (avg, list(y)) #dictionary with gene as key, tuple of avd std dev and list of means as value

    x = [i for i in range(1, groupby.shape[0] + 1)]
    fig = plt.figure(figsize=(20, 2))
    ax = fig.add_subplot(111)
    ax.set_xticks(x)
    ax.errorbar(x, y, yerr=y_err, fmt='o')
    ax.set_xlabel('CpG #', fontsize=12)
    ax.set_ylabel('Average Methylation Ratio', fontsize=12)
    fig.savefig(cancer_type + row[3] + '_cds_CpGs.png')





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

