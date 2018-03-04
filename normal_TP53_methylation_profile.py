import numpy as np
import matplotlib.pyplot as plt
import pybedtools
import pandas as pd
import csv
import fileinput


def plot_profile(file):
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

    x = [i for i in range(1, groupby.shape[0] + 1)]
    fig = plt.figure(figsize=(20, 2))
    ax = fig.add_subplot(111)
    ax.set_xticks(x)
    ax.errorbar(x, y, yerr=y_err, fmt='o')
    ax.set_xlabel('CpG #', fontsize=12)
    ax.set_ylabel('Average Methylation Ratio', fontsize=12)
    fig.savefig('TP53_normal_methylation_profile.png')

# a = pybedtools.BedTool('all_normal_tissues_WGBS_cds.bed')
# a.intersect('TP53_cds.bed').saveas('TP53_normal_cds_CpGs.bed')

# #add headers
# headers = 'chromosome start stop id methylation_ratio'.split()
#
# for line in fileinput.input(['TP53_normal_cds_CpGs.bed'], inplace=True):
#     print ('\t'.join(headers))
#     print (line)
plot_profile('TP53_normal_cds_CpGs.bed')