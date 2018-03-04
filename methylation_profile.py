from multiprocessing import Pool
import pybedtools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import pickle

#script to generate methylation profile of coding regions of a given gene
#profile shows for each CpG in coding region of the gene the methylation ratio observed in normal tissue, the average 
#methylation ratio observed across many cancer samples along with std dev, and sites where a mutation occurred if it was seen in at
#least one cancer sample, along with the frequency of samples it was observed in

#INPUTS:
#File 1 --> Whole genome bisulfite sequencing from cancer tissue(many samples)
#File 2 --> Whole genome bisulfite sequencing from matched normal tissue(one sample)
#File 3 --> Simple somatic mutation data from same cancer type(many samples)



#input WGBS file from cancer tissue
file = '/Users/khandekara2/Documents/methylationProject/01_data/MALY_prevBase_TP53cds.tsv'
df = pd.read_csv(file, sep='\t')
locations = []
for chrom, start, stop in zip(list(df['chromosome']), list(df['start']), list(df['stop'])):
    locations.append((str(chrom), int(start), int(stop)))
df['location'] = locations
grouped = df.groupby('location')
groupby = grouped.agg({'methylation_ratio': [np.mean, np.std], 'prev_base': 'max'})
groupby.reset_index(inplace=True)
groupby.columns = groupby.columns.droplevel(0)

y = groupby['mean'] #mean ratio from cancer
y_err = groupby['std'] #standard deviation from cancer samples
x = np.array([i for i in range(1, groupby.shape[0] + 1)]) # ticks for each CpG in coding region
x_labels = list(set([i[1] for i in df.location])) #coordinate of cytosine as labels
# x_labels = [i + 'CG' for i in groupby['max']]
x_labels.sort()

#input WGBS file from normal tissue
df2 = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/TP53_normal_cds_CpGs.tsv', sep='\t')
c = list(df2.start)
c = [i for i in c if i in x_labels]
c.sort() # now normal and cancer coordinates are matching
y2 = [] #methylation ratio from normal tissue 
for e in c:
    r = df2['methylation_ratio'].where(df2['start'] == e)
    # r is a series with all Nan's except for the value we want
    for a in r:
        if not np.isnan(a):
            y2.append(a)

#add mutated sites to profile
mutation_dict = pickle.load(open('/Users/khandekara2/Documents/methylationProject/01_data/mutation_dict_pickles/MALY_mutation_dict.pickle', 'rb'))
mutations = [] #coordinate of mutated site
mut_ratios = [] #metthylation ratio of mutated site if available
frequencies = [] #frequency of mutation(# of samples it occurred in)
for loc, sample, ratio in zip(df.location, df.id, df.methylation_ratio):
#     loc = (str(chrom), int(start), int(stop))
    if loc in mutation_dict and loc[1] not in mutations:
        if sample in mutation_dict[loc]:
            mut_ratios.append(ratio)
        else:
            mut_ratios.append(1.1)
        
        mutations.append(loc[1])
        frequencies.append(len(mutation_dict[loc]))
               
# mut_ratios = [0.0 for _ in mutations]
print (np.array(list(mutations)))
print (mut_ratios)

#construct plot/profile
%matplotlib inline
fig = plt.figure(figsize=(20, 2))
ax = fig.add_subplot(111)
ax.set_xticks(x)
ax.set_xticklabels(x_labels, rotation='vertical', fontsize=12)
ax.set_ylim((0.0, 1.1))
ax.errorbar(x, y, yerr=y_err, fmt='o')
ax.set_xlabel('Coordinate', fontsize=14)
ax.set_ylabel('Average Methylation Ratio', fontsize=14)
ax.scatter(x, y2, color='red')
# ax.scatter(np.array(list(mutations)), np.array(mut_ratios), color='green')
trans = ax.get_xaxis_transform()
for u, v, f in zip(mutations, mut_ratios, frequencies):
    ax.annotate(str(f), xy=(x_labels.index(u) + 1, v), xycoords=trans, arrowprops=dict(facecolor='black', shrink=0.5,))
ttl = ax.title
ttl.set_position([.5, 1.2])
ax.set_title('Methylation Profile of TP53 Coding Region in MALY', fontsize=16)
# plt.gcf().subplots_adjust(bottom=0.15)
fig.savefig('/Users/khandekara2/Documents/methylationProject/03_results/MALY_TP53_methylation_profile.png', bbox_inches='tight')