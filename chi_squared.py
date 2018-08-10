import pandas as pd
from scipy.stats import chi2_contingency
import sys
import decimal
import csv
import numpy as np

with open('chi2_stats.tsv', 'w') as csvout:
    writer = csv.writer(csvout, delimiter='\t')
    effect_sizes = []
    df = pd.read_csv('/data/khandekara2/imputation/mutations_reconstructed.tsv', sep='\t')
    prefix = 'MALY_'
    suffix = '_WGBS.bed.single.non_mutated'
    for sample in df['id'].unique():
        #get counts for mutated sites
        mut = df[df['id'] == sample]
        mm = mut[mut['prior_ratio'] >= 0.5].shape[0]
        um = mut[mut['prior_ratio'] < 0.5].shape[0]
        #get counts for non-mutated sites
        non_mut = pd.read_csv(prefix + sample + suffix, sep='\t')
        nm = non_mut[non_mut.iloc[:, 3] >= 0.5].shape[0]
        nn = non_mut[non_mut.iloc[:, 3] < 0.5].shape[0]
        obs = np.array([[mm, nm], [um, nn]])
        chi2, p, dof, expected = chi2_contingency(obs)
        effect_size = np.sqrt(chi2 / (mm + nm  + um + nn)) #Cramers V
        effect_sizes.append(effect_size)
        writer.writerow([sample, round(p, 2), round(effect_size, 2)])
