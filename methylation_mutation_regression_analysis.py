
from __future__ import division
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys

mutToFrequency = {} #key is bin, value is frequency of mutated CpG sites in this bin
totalCpGFrequency = {} #key is bin, value is total frequency of CpG sites in this bin

#takes a list of methylation ratios and "bins" them while keeping track of the frequency in each bin.
#ratios of exactly 0 and exactly 1 are separate bins, and the rest are in increments of 0.1 from 0 to 1 for a total of 12 bins
#for example, a bin of 0.25 represents all methylation ratios falling between 0.2 and 0.3

def bin(binToFrequency, ratios):
    for ratio in ratios:
        if ratio == 0: #0
            if 0 not in binToFrequency.keys():
                binToFrequency[0] = 0
            binToFrequency[0] += 1

        elif ratio > 0 and ratio <= 0.1: #0.05
            if 0.05 not in binToFrequency.keys():
                binToFrequency[0.05] = 0
            binToFrequency[0.05] += 1

        elif ratio > 0.1 and ratio <= 0.2: #0.15
            if 0.15 not in binToFrequency.keys():
                binToFrequency[0.15] = 0
            binToFrequency[0.15] += 1

        elif ratio > 0.2 and ratio <= 0.3: #0.25
            if 0.25 not in binToFrequency.keys():
                binToFrequency[0.25] = 0
            binToFrequency[0.25] += 1

        elif ratio > 0.3 and ratio <= 0.4: #0.35
            if 0.35 not in binToFrequency.keys():
                binToFrequency[0.35] = 0
            binToFrequency[0.35] += 1

        elif ratio > 0.4 and ratio <= 0.5: #0.45
            if 0.45 not in binToFrequency.keys():
                binToFrequency[0.45] = 0
            binToFrequency[0.45] += 1

        elif ratio > 0.5 and ratio <= 0.6: #0.55
            if 0.55 not in binToFrequency.keys():
                binToFrequency[0.55] = 0
            binToFrequency[0.55] += 1

        elif ratio > 0.6 and ratio <= 0.7: #0.65
            if 0.65 not in binToFrequency.keys():
                binToFrequency[0.65] = 0
            binToFrequency[0.65] += 1

        elif ratio > 0.7 and ratio <= 0.8: #0.75
            if 0.75 not in binToFrequency.keys():
                binToFrequency[0.75] = 0
            binToFrequency[0.75] += 1

        elif ratio > 0.8 and ratio <= 0.9: #0.85
            if 0.85 not in binToFrequency.keys():
                binToFrequency[0.85] = 0
            binToFrequency[0.85] += 1

        elif ratio > 0.9 and ratio < 1: #0.95
            if 0.95 not in binToFrequency.keys():
                binToFrequency[0.95] = 0
            binToFrequency[0.95] += 1

        else: #ratio == 1
            if 1 not in binToFrequency.keys():
                binToFrequency[1] = 0
            binToFrequency[1] += 1

def plot():
    xs = np.array(list(mutToFrequency.keys()))
    ys = np.array([mutToFrequency[x] / totalCpGFrequency[x] * 100000  for x in mutToFrequency.keys()])  #frequency of mutations per megabase
    slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)
    line = slope * xs + intercept
    plt.xticks([0, 0.25, 0.5, 0.75, 1])
    plt.plot(xs, ys, 'o', xs, line, '-')
    plt.ylabel('Mutations per MB of CpG dinucleotides')
    plt.xlabel('Fraction CpGs Methylated')
    plt.gcf().text(1.0, 1.0, "Slope: %.9f" % slope, fontsize=12)
    plt.gcf().text(1.0, 0.9, "P-value: %.9f" % p_value, fontsize=12)
    plt.gcf().text(1.0, 0.8, "R-squared: %.2f" % r_value ** 2, fontsize=12)
    if file1.endswith('non_cds'):
        plt.title(file2.split('-')[0] + ' Non-Coding Region')
        plt.savefig(file2.split('-')[0] + '_Non-CodingRegion.png', bbox_inches='tight')
        plt.clf()
    else:
        plt.title(file2.split('-')[0] + ' Coding Region')
        plt.savefig(file2.split('-')[0] + '_CodingRegion.png', bbox_inches='tight')
        plt.clf()


#Input file 1: bed file containing all CpG's assayed and their mean methylation ratio from normal tissue
#Input file 2: bed file containing all CpG sites where mutations occurred


# all_sites = ['GSM1010984_UCSD.Gastric.Bisulfite-Seq.STL003.wig.bed.non_cds', 'GSM1010984_UCSD.Gastric.Bisulfite-Seq.STL003.wig.bed.cds', 'normal_colon_WGBS_processed_expanded.bed.non_cds', 'normal_colon_WGBS_processed_expanded.bed.cds', 'GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig.bed.non_cds', 'GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig.bed.cds', 'GSM983647_UCSD.Lung.Bisulfite-Seq.STL002.wig.bed.non_cds', 'GSM983647_UCSD.Lung.Bisulfite-Seq.STL002.wig.bed.cds', 'GSM916052_BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig.bed.non_cds', 'GSM916052_BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig.bed.cds', 'GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.bed.non_cds', 'GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.bed.cds']
# mutated_sites = ['GASTRIC-COMBINEDmethylation_mutation.bed', 'GASTRIC-COMBINEDmethylation_mutation_cds.bed', 'COAD-USmethylation_mutation.bed', 'COAD-USmethylation_mutation_cds.bed', 'LIVER-COMBINEDmethylation_mutation.bed', 'LIVER-COMBINEDmethylation_mutation_cds.bed', 'LUNG-COMBINEDmethylation_mutation.bed', 'LUNG-COMBINEDmethylation_mutation_cds.bed', 'CLLE-ESmethylation_mutation.bed', 'CLLE-ESmethylation_mutation_cds.bed', 'OV-USmethylation_mutation.bed', 'OV-USmethylation_mutation_cds.bed']
all_sites = ['GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.wig.bed.non_cds', 'GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.wig.bed.cds', 'GSM983651_UCSD.Pancreas.Bisulfite-Seq.STL003.wig.bed.non_cds', 'GSM983651_UCSD.Pancreas.Bisulfite-Seq.STL003.wig.bed.cds', 'GSM1010979_UCSD.Thymus.Bisulfite-Seq.STL001.wig.bed.non_cds', 'GSM1010979_UCSD.Thymus.Bisulfite-Seq.STL001.wig.bed.cds']
mutated_sites = ['BRCA-USmethylation_mutation.bed', 'BRCA-USmethylation_mutation_cds.bed', 'PANCREAS-COMBINEDmethylation_mutation.bed', 'PANCREAS-COMBINEDmethylation_mutation_cds.bed', 'MALY-DEmethylation_mutation.bed', 'MALY-DEmethylation_mutation_cds.bed']

for file1, file2 in zip(all_sites, mutated_sites):
    df = pd.read_csv(file1, sep='\t', header=None)
    mut = pd.read_csv(file2, sep='\t', header=None)
    num_samples = mut.iloc[:, 3].nunique()
    means = df.iloc[:, 4] #normal methylation ratios
    means2 = mut.iloc[:, 4]

    bin(totalCpGFrequency, means)
    bin(mutToFrequency, means2)
    plot()