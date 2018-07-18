#simple linear regression to model the relationship between methylation ratio and mutation frequency in various cancers
#goal is to determine how much methylation expalins the variation in coding and non-coding regions compare the models derived from coding regions vs. all other regions
#normal WGBS data from various tissues is used in conjunction with mutation data of cancers originating in same tissue
#version which plots graphs on one figure

from bin_ratios import *
from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np
from scipy import stats
import sys
import decimal
%matplotlib inline
import pybedtools
import os

# matching of WGBS data from normal tissue to pooled mutation in the same tissue
wgbs_files = ['GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig.processed',  'GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.wig.bed.processed', 'GSM983647_UCSD.Lung.Bisulfite-Seq.STL002.wig.bed.processed', 'GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.bed.processed',  'GSM983651_UCSD.Pancreas.Bisulfite-Seq.STL003.wig.bed.processed', 'GSM1010984_UCSD.Gastric.Bisulfite-Seq.STL003.wig.bed.processed', 'GSM1010979_UCSD.Thymus.Bisulfite-Seq.STL001.wig.bed.processed', 'normal_colon_WGBS_processed_expanded.bed.processed', 'GSM916052_BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig.bed.processed']
mutation_files = ['LIVER-COMBINED_mutation.bed.all.noDuplicates', 'BRCA-US_mutation.bed.all.noDuplicates', 'LUNG-COMBINED_mutation.bed.all.noDuplicates', 'OV-US_mutation.bed.all.noDuplicates', 'PANCREAS-COMBINED_mutation.bed.all.noDuplicates', 'GASTRIC-COMBINED_mutation.bed.all.noDuplicates', 'MALY-DE_mutation.bed.all.noDuplicates', 'COAD-US_mutation.bed.all.noDuplicates', 'CLLE-ES_mutation.bed.all.noDuplicates', ]

#create input files for regression plots
for meth, mut in zip(wgbs_files, mutation_files:
	a = pybedtools.BedTool(meth)
	a.intersect(mut, wa=True wb=True).saveas(mut.split('_')[0] + 'methylation_mutation.bed') #non-cds

mutated_sites = ['breast_methylation_mutation.bed', 'pancreas_methylation_mutation.bed',  'thymus_methylation_mutation.bed', 'lung_methylation_mutation.bed',  'ovary_methylation_mutation.bed', 'CLLE_methylation_mutation.bed', 'esophagus_methylation_mutation.bed', 'gastric_methylation_mutation.bed', 'colon_methylation_mutation.bed']
#construct regression plots using our method
if __name__ == "__main__":
	fig = plt.figure(figsize=(3,6))
	for i, (file1, file2) in enumerate(zip(wgbs_files, mutated_sites)):
        df = pd.read_csv(file1, sep='\t', header=None) #file containing ratios from normal tissue
        mut = pd.read_csv(file2, sep='\t', header=None) #file containing C>T mutations in CpG's along with their corresponding rato in normal tissue
		num_samples = mut.iloc[:, 3].nunique()
		xToY = {0:[], 0.05:[], 0.15:[], 0.25:[], 0.35:[], 0.45:[], 0.55:[], 0.65:[], 0.75:[], 0.85:[], 0.95:[], 1:[]} #key is bin, value is list of length num_samples
		totalCpGFrequency = {} #key is bin, value is total frequency of CpG sites in this bin
		bin(totalCpGFrequency, ratios)
		for sample in mut.iloc[:, 3].unique()
			mut = mut[mut['id'] == sample]
        	ratios = df.iloc[:, 4]
        	mut_ratios = mut.iloc[:, 4]
			mutToFrequency = {} #key is bin, value is frequency of mutated CpG sites in this bin
        	bin(mutToFrequency, mut_ratios)
			for key, value in mutToFrequency.items():
				xToY[key].append(value / (totalCpGFrequency[key] / 1000000) #normalized count of mutations for this specific sample is appended to appropiate bin
		xs = np.array(sorted(list(mutToFrequency.keys()))) * num_samples
		ys = []
		for key in sorted(xToY.keys()):
			ys += xToY[key]
		print (len(xs) == len(ys))






#create matrix
for f in wgbs_files:
	for m in mutation_files:

# print(len(methylation_files) == len(mutation_files) == len(cds_files))
