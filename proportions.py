#script to calculate what proportion of methylated CpG's fall into each bin
from __future__ import division
import pandas as pd
import sys
import csv


filename = sys.argv[1]
df = pd.read_csv(filename,  sep='\t')

with open('proportion_results.csv', 'w') as f:
    writer = csv.writer(f)
    #write header
    writer.writerow(['sample', 'lowly', 'low-intermediately', 'intermediately', 'high-intermediately', 'highly methylated'])
    biseq_samples = df['id'].unique()
    lowly = 0
    low_intermediately = 0
    intermediately = 0
    high_intermediately = 0
    highly_methylated = 0
    for sample in biseq_samples:
        sub_biseq = df[df['id'] == sample]
        total_sites = sub_biseq.shape[0]

        for ratio in sub_biseq['methylation_ratio']:
            if ratio > 0 and ratio <= 0.2:
                lowly +=1
            if ratio > 0.2 and ratio <= 0.4:
                low_intermediately += 1
            if ratio > 0.4 and ratio <= 0.6:
                intermediately += 1
            if ratio > 0.6 and ratio <= 0.8:
                high_intermediately += 1
            if ratio > 0.8 and ratio <= 1:
                highly_methylated += 1

        #transform into proportions
        lowly = (lowly / total_sites) * 100
        low_intermediately = (low_intermediately / total_sites) * 100
        intermediately = (intermediately / total_sites) * 100
        high_intermediately = (high_intermediately / total_sites) * 100
        highly_methylated = (highly_methylated / total_sites) * 100

        writer.writerow([sample, lowly, low_intermediately, intermediately, high_intermediately, highly_methylated])


















