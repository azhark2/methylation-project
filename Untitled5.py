
# coding: utf-8

# In[135]:

import pickle as cPickle
import numpy

WGBS_avgBeta = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/WGBS_avgBeta.pickle", "rb"))
# # cancer_cds_overlap = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/cancer_cds_overlap", "rb"))
# cds_overlap = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/cds_overlap", "rb"))


# In[60]:

print WGBS_avgBeta.keys()



# In[61]:

print WGBS_avgBeta['CLLE'].items()[-5]


# In[62]:

print len(cancer_cds_overlap)


# In[63]:

print len(cds_overlap)


# In[64]:

ratios = []
for chrom, location in cancer_cds_overlap:
    sum = 0
    count = 0
    if (chrom, location) in WGBS_avgBeta['MALY']:
        sum += WGBS_avgBeta['MALY'][(chrom, location)]
        count += 1
    if (chrom, location) in WGBS_avgBeta['PBCA']:
        sum += WGBS_avgBeta['PBCA'][(chrom, location)]
        count += 1
    if (chrom, location) in WGBS_avgBeta['CLLE']:
        sum += WGBS_avgBeta['CLLE'][(chrom, location)]
        count += 1
    avg = sum / 3
    ratios.append(avg)
print ratios
    


# In[65]:

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
plt.hist(ratios, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))


# In[139]:

cancer_overlaps = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/cancer_overlaps.pickle", "rb"))


# In[140]:

print cancer_overlaps.items()


# In[144]:

values = []
for location in cancer_overlaps['PBCA'][0]:
    values.append(WGBS_avgBeta['PBCA'][location])
plt.hist(values, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Methylation Ratio')
plt.ylabel('Frequency of Overlaps')
plt.title('PBCA')


# In[142]:

print len(cancer_overlaps['PBCA'][0])


# In[143]:

print len(cancer_overlaps['CLLE'][0])


# In[147]:

values = []
for location in cancer_overlaps['MALY'][0]:
    if location in WGBS_avgBeta['MALY']:
        values.append(WGBS_avgBeta['MALY'][location])
plt.hist(values, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Methylation Ratio')
plt.ylabel('Frequency of Overlaps')
plt.title('MALY')


# In[148]:

print len(values)


# In[72]:

values = []
for location in cancer_overlaps['CLLE'][0]:
    if location in WGBS_avgBeta['CLLE']:
        values.append(WGBS_avgBeta['CLLE'][location])
plt.hist(values, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))


# In[75]:

import pandas as pd
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
biseq = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_CLLE_cds.tsv', sep='\t')
biseq.hist(column='methylation_ratio', bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))


# In[77]:

values = []
for location in cancer_overlaps['CLLE'][0]:
    if location in WGBS_avgBeta['CLLE']:
        print (WGBS_avgBeta['CLLE'][location])


# In[128]:

df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/03_results/results.csv')


# In[83]:

#GSEA prep
f = open('gene_ids', 'w')
mutation = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv', sep='\t')
mutation = mutation.dropna(subset=['gene_affected'])

locations = []
for chrom, coord in zip(list(mutation['chromosome']), list(mutation['chromosome_start'])):
    locations.append(('chr' + str(chrom), str(coord)))

pbca_overlaps = set(cancer_overlaps['PBCA'][0])
mutation['location'] = locations
gene_ids = set([])
for overlap, gene in zip(mutation['location'], mutation['gene_affected']):
    if overlap in pbca_overlaps:
        gene_ids.add(gene)

for gene in gene_ids:   
    f.write(gene)
    f.write('\n')
        
        
        
    
    
    
    

df = df[df['cancer_type'].str.contains('PBCA')]


# In[123]:

df.head()


# In[129]:

df = df[df['cancer_type'].str.contains('MALY')]


# In[130]:

import matplotlib.pyplot as plt
import pylab

y = df['num_cds_cpg_sites_mutated']
x = df['num_cds_sites_methylated']

plt.scatter(x,y)
# plt.title('Number of CpG Sites Mutated vs Number of Sites Methylated in PBCA')
plt.ylabel('Number of CpG Sites Mutated in Protein Coding Regions')
plt.xlabel('Number of CpG Sites Methylated in Protein Coding Regions')

plt.show()


# In[154]:

import collections
counter=collections.Counter(cancer_overlaps['MALY'][0])
print(counter)


# In[155]:

print(counter.most_common(5))


# In[157]:

print len(WGBS_avgBeta['PBCA'].keys())


# In[158]:

a = ('1', '1000')
b = ('1', '1000')
print (a == b)


# In[165]:

a = set([1, 2, 3])
b = set([3, 5, 6])
c = a and b 
type(c)


# In[166]:

c = set.intersection(a, b)
print c


# In[ ]:



