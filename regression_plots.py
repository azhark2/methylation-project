from bin_ratios import *
from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
import numpy as np
from scipy import stats
import sys
import decimal

wgbs = ['GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig.bed',  'GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.wig.bed', 'GSM983647_UCSD.Lung.Bisulfite-Seq.STL002.wig.bed', 'GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.bed',  'GSM983651_UCSD.Pancreas.Bisulfite-Seq.STL003.wig.bed', 'GSM983649_UCSD.Esophagus.Bisulfite-Seq.STL003.wig.bed' 'GSM1010984_UCSD.Gastric.Bisulfite-Seq.STL003.wig.bed', 'GSM1010979_UCSD.Thymus.Bisulfite-Seq.STL001.wig.bed',  'GSM916052_BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig.bed']
mutations = ['Liver_methylation_mutation.bed', 'Breast_methylation_mutation.bed', 'Lung_methylation_mutation.bed', 'Ovary_methylation_mutation.bed', 'Pancreas_methylation_mutation.bed', 'Esophagus_methylation_mutation.bed', 'Stomach_methylation_mutation.bed', 'MALY_methylation_mutation.bed', 'CLLE_methylation_mutation.bed']
fig = plt.figure(figsize=(6,4))
fig2 = plt.figure(figsize=(6,4))
for i, (normal, mutated) in enumerate(zip(wgbs, mutations)):
    df = pd.read_csv(normal, sep='\t')
    mut = pd.read_csv(mutated, sep='\t')
    title = mutated.split('_')[0]
    num_samples = mut['id'].nunique()
    xToY = {0:[], 0.05:[], 0.15:[], 0.25:[], 0.35:[], 0.45:[], 0.55:[], 0.65:[], 0.75:[], 0.85:[], 0.95:[], 1:[]} #key is bin, value is list of length num_samples
    totalCpGFrequency = {} #key is bin, value is total frequency of CpG sites in this bin
    ratios = df.iloc[:, 4]
    bin(totalCpGFrequency, ratios)
    num_mutations = []
    correlations = []
    num_samples = mut['id'].nunique()
    for sample in mut['id'].unique():
        sub_mut = mut[mut['id'] == sample] #go through each sample one by one
        num_mut = sub_mut.shape[0]
        num_mutations.append(num_mut)
        sub_mut_ratios = sub_mut['methylation_ratio']
        mutToFrequency = {} #key is bin, value is frequency of mutated CpG sites in this bin
        bin(mutToFrequency, sub_mut_ratios)
        x_val = []
        y_val = []
        for key in sorted(mutToFrequency.keys()):
            xToY[key].append((mutToFrequency[key] / num_mut) / (totalCpGFrequency[key] / 1000000)) #normalized count of mutations for this specific sample is appended to appropiate bin
            x_val.append(key)
            y_val.append((mutToFrequency[key] / num_mut) * 100 / (totalCpGFrequency[key] / 1000000))
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_val, y_val)
        correlations.append(r_value ** 2)

    #plot histogram of correlations for each sample (we have many samples for each cancer type)
    ax2 = fig2.add_subplot(3, 3, i+1)
    ax2.hist(correlations, edgecolor="k")
    ax2.set_title(title)
    if i+1 == 1 or i+1 == 4 or i+1 == 7:
        ax2.set_ylabel('Frequency of\n Samples', fontsize=7)
    ax2.set_xlabel(r'$R^2$')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    xs = []
    ys = []
    for key in sorted(xToY.keys()):
        for k in xToY[key]:
            if k != 0:
                xs.append(key)
                ys.append(k)

    print (len(xs) == len(ys))
    xs = np.array(xs)
    ys = np.array(ys)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)
#     print(slope, intercept, r_value ** 2, p_value, std_err)
    line = slope * xs + intercept

    #construct plot
    ax = fig.add_subplot(3, 3, i+1)
    if i+1 == 1 or i+1 == 4 or i+1 == 7:
        ax.set_xlabel('Fraction CpGs Methylated', fontsize=7)
        ax.set_ylabel('% CpG Mutations\nper MB', fontsize=7, multialignment='center')
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.plot(xs, ys, 'o', xs, line, '-', markersize=2)
    ax.set_title(title)
    #ax.text(0.2, 0.85, "P-value: %.2E" % decimal.Decimal(p_value), fontsize=6, verticalalignment='center', horizontalalignment='left', clip_on=True, transform=ax.transAxes)
    ax.text(0.1, 0.75, r'$R^2$: %.2f' % r_value ** 2, fontsize=7, verticalalignment='center', horizontalalignment='center', clip_on=True, transform=ax.transAxes)
    #ax.text(0.2, 0.75, "Slope: %.2E" % decimal.Decimal(slope), fontsize=6, verticalalignment='center', horizontalalignment='left', clip_on=True, transform=ax.transAxes)
    ax.text(0.9, 0.75, "n = " + str(num_samples), fontsize=7, verticalalignment='center', horizontalalignment='left', clip_on=True, transform=ax.transAxes)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.tight_layout()
fig2.tight_layout()
fig.savefig('scatterplots.png', dpi=600, bbox_inches="tight")
fig2.savefig('correlation_distributions.png', dpi=600, bbox_inches="tight")
