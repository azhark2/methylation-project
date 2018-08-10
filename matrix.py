from __future__ import division
from bin_ratios import *
import pandas as pd
import numpy as np
from scipy import stats
import sys
import decimal
import csv
import pybedtools
import os

wgbs_files = ['GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig.bed',  'GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.wig.bed', 'GSM983647_UCSD.Lung.Bisulfite-Seq.STL002.wig.bed', 'GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.bed',  'GSM983651_UCSD.Pancreas.Bisulfite-Seq.STL003.wig.bed', 'GSM983649_UCSD.Esophagus.Bisulfite-Seq.STL003.wig.bed', 'GSM1010984_UCSD.Gastric.Bisulfite-Seq.STL003.wig.bed', 'GSM1010979_UCSD.Thymus.Bisulfite-Seq.STL001.wig.bed',  'GSM916052_BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig.bed']
mutations = ['Liver_methylation_mutation.bed', 'Breast_methylation_mutation.bed', 'Lung_methylation_mutation.bed', 'Ovary_methylation_mutation.bed', 'Pancreas_methylation_mutation.bed', 'Esophagus_methylation_mutation.bed', 'Stomach_methylation_mutation.bed', 'MALY_methylation_mutation.bed', 'CLLE_methylation_mutation.bed']

with open('matrix.tsv', 'w') as csvout:
    writer = csv.writer(csvout, delimiter='\t')
    for normal in wgbs_files:
        correlations = []
        a = pybedtools.BedTool(normal)
        for mutated in mutations:
            a.intersect(mutated, wa=True, wb=True).saveas(normal.split('.')[0] + '_' + mutated.split('_')[0] + '.matched')
            df = pd.read_csv(normal, sep='\t')
            mut = pd.read_csv(normal.split('.')[0] + '_' + mutated.split('_')[0] + '.matched', sep='\t')
            xToY = {0:[], 0.05:[], 0.15:[], 0.25:[], 0.35:[], 0.45:[], 0.55:[], 0.65:[], 0.75:[], 0.85:[], 0.95:[], 1:[]} #key is bin, value is list of length num_samples
            totalCpGFrequency = {} #key is bin, value is total frequency of CpG sites in this bin
            ratios = df.iloc[:, 4]
            bin(totalCpGFrequency, ratios)
            num_mutations = []
            num_samples = mut.iloc[:, 9].nunique()
            max_correlation = 0
            for sample in mut.iloc[:, 9].unique():
                sub_mut = mut[mut.iloc[:, 9] == sample] #go through each sample one by one
                num_mut = sub_mut.shape[0]
                if num_mut >= 50:
                    num_mutations.append(num_mut)
                    sub_mut_ratios = sub_mut.iloc[:, 4]
                    mutToFrequency = {} #key is bin, value is frequency of mutated CpG sites in this bin
                    bin(mutToFrequency, sub_mut_ratios)
                    x_val = []
                    y_val = []
                    for key in sorted(mutToFrequency.keys()):
                        xToY[key].append((mutToFrequency[key] / num_mut) / (totalCpGFrequency[key] / 1000000)) #normalized count of mutations for this specific sample is appended to appropiate bin
                        x_val.append(key)
                        y_val.append((mutToFrequency[key] / num_mut) * 100 / (totalCpGFrequency[key] / 1000000))
                    slope, intercept, r_value, p_value, std_err = stats.linregress(x_val, y_val)
                    if r_value ** 2 > max_correlation:
                        max_correlation = r_value ** 2
            os.remove(normal.split('.')[0] + '_' + mutated.split('_')[0] + '.matched')
            correlations.append(r_value ** 2)
        writer.writerow(correlations)
