
# coding: utf-8

# In[1]:

import _pickle as cPickle


# In[2]:


sitecount = cPickle.load(open('siteCount.pickle', 'rb'))


# In[3]:

print (sitecount.items())


# In[5]:

all_betaValues = cPickle.load(open('all_beta_values .pickle', 'rb'))


# In[6]:

gene_dict = cPickle.load(open('/Users/khandekara2/Documents/methylationProject/02_code/gene_dict.pickle ', 'rb'))


# In[7]:

gene_dict2 = cPickle.load(open('/Users/khandekara2/Documents/methylationProject/02_code/gene_dict2.pickle ', 'rb'))


# In[8]:

sample_count = cPickle.load(open('/Users/khandekara2/Documents/methylationProject/02_code/sample_count.pickle', 'rb'))
print (len(sample_count))
files = ['https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BLCA-US/meth_array.BLCA-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-US/meth_array.BRCA-US.tsv.gz',
       'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/CESC-US/meth_array.CESC-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/CLLE-ES/meth_array.CLLE-ES.tsv.gz',
    'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LAML-US/meth_array.LAML-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PAEN-AU/meth_array.PAEN-AU.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/GBM-US/meth_array.GBM-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/SKCM-US/meth_array.SKCM-US.tsv.gz',
             'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/UCEC-US/meth_array.UCEC-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/STAD-US/meth_array.STAD-US.tsv.gz',
             'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-US/meth_array.PRAD-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-CA/meth_array.PRAD-CA.tsv.gz',
             'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PBCA-DE/meth_array.PBCA-DE.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/OV-US/meth_array.OV-US.tsv.gz',
             'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/OV-AU/meth_array.OV-AU.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/HNSC-US/meth_array.HNSC-US.tsv.gz',
             'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/READ-US/meth_array.READ-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/THCA-US/meth_array.THCA-US.tsv.gz',
             'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LGG-US/meth_array.LGG-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LIHC-US/meth_array.LIHC-US.tsv.gz',
         'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/GBM-US/meth_array.GBM-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LAML-US/meth_array.LAML-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/KIRP-US/meth_array.KIRP-US.tsv.gz',
         'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/KIRC-US/meth_array.KIRC-US.tsv.gz']
print (len(files))


# In[9]:

# sample_dict = {type: 0 for type in sitecount.keys()} #dictionary that maps cancer type to # of samples
# for count in sample_count:
    
sample_dict = {}
# print (sample_dict)    
for f, samples in zip(files, sample_count):
    if f[len(f) - 14:-10] in sample_dict:
        sample_dict[f[len(f) - 14:-10]] += samples
    else:
        sample_dict[f[len(f) - 14:-10]] = samples

 
print (len(sample_dict.keys()))
print (sample_dict.items())
    


# In[10]:

plt.bar(range(len(sample_dict)), sample_dict.values(), align='center')
plt.xticks(range(len(sample_dict)), sample_dict.keys(), rotation = 'vertical')
plt.yticks(np.arange(0, 1300, 100))
plt.xlabel("Type of Cancer")
plt.ylabel("# of samples")

plt.show()


# In[11]:

sorted_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse = True) #list of tuples of genes and the proportion of methylated sites(from highest to lowest)


# In[12]:

print (sorted_genes[:20])


# In[13]:

import pandas as pd
df1 = pd.read_csv('Census_allThu Feb  9 18_47_50 2017.csv', header = 0) #Cosmic cancer genes
cosmicGenes = set(df1.iloc[:, 0])
# print (cosmicGenes[:20])


# In[14]:

topCosmicGenes = []
topCosmicCounts = []
sorted_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse = True) #list of tuples of genes and the proportion of methylated sites(from highest to lowest)
genes = [pair[0] for pair in sorted_genes]
proportion_methylated = [pair[1] for pair in sorted_genes]


for gene, count in zip(genes, proportion_methylated):
    if gene in cosmicGenes:
        topCosmicGenes.append(gene)
        topCosmicCounts.append(count)


# In[15]:

top = topCosmicGenes[:20]
print (topCosmicGenes[:20])
print (topCosmicCounts[:20])


# In[16]:

totals = []
for gene in top:
    totals.append(gene_dict2[gene])
    


# In[17]:

print (totals)


# In[18]:

print (topCosmicGenes[:20])


# In[19]:

for key in sitecount.keys():
    sitecount[key][0] = int(sitecount[key][0] / sample_dict[key])
    sitecount[key][1] = int(sitecount[key][1] / sample_dict[key])
print(sitecount.items())

sorted_count = sorted(sitecount.items(), key = lambda x: x[1][1], reverse = True)
print(sorted_count)


# In[20]:

top_cancers = ['LIHC', 'CLLE', 'PRAD', 'PBCA', 'UCEC']
top_exon_sites = [2841, 4491, 3867, 4040, 2989]
top_total_sites = [107589, 165721, 138736, 155467, 93422]


# In[21]:

cancer_types = []
exon_sites = []
total_sites = []
for key, value in sitecount.items():
    cancer_types.append(key)
    exon_sites.append(value[1])
    total_sites.append(value[0])

cancer_types[6] = 'GBM'
cancer_types[7] = 'LGG'
cancer_types[11] = 'OV'
print (cancer_types)
print (exon_sites[:5])
print (total_sites[:5])
print (len(cancer_types))
    
    
    


# In[ ]:



    


# In[22]:

import matplotlib.pyplot as plt
# plt.rcParams['figure.figsize'] = (20.0, 10.0)
# import matplotlib.ticker as ticker
import numpy as np
get_ipython().magic('matplotlib inline')
ind = np.arange(5) #np.arange(len(cancer_types))
width = 0.35
# figure = pylab.figure()
# ax = figure.add_subplot(1,1,1)
fig, ax = plt.subplots(nrows = 1, ncols=1)
rects1 = ax.bar(ind, top_exon_sites, color = 'r')
rects2 = ax.bar(ind + width, top_total_sites, width, color='y')
ax.set_xticks(ind + width/2)


# tick_spacing = 20
# ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
ax.set_xticklabels(top_cancers, rotation = 'vertical')
plt.tight_layout()
# def autolabel(rects):
#     """
#     Attach a text label above each bar displaying its height
#     """
#     for rect in rects:
#         height = rect.get_height()
#         ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
#                 '%d' % int(height),
#                 ha='center', va='bottom')

# autolabel(rects1)
# autolabel(rects2)


# In[23]:

import matplotlib.pyplot as plt

import matplotlib.ticker as ticker
import numpy as np
get_ipython().magic('matplotlib inline')
ind = np.arange(len(cancer_types))
# plt.figure(figsize=(20,10))
# fig.add_axes()
width = 0.35
fig, ax = plt.subplots(nrows = 1, ncols=1)

rects1 = ax.bar(np.arange(len(cancer_types)), exon_sites, color = 'r')
# rects2 = ax.bar(ind + width, total_sites, width, color='y')
ax.set_xticks(ind)

# tick_spacing = 20
# ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
ax.set_xticklabels(cancer_types, rotation = 'vertical')
plt.title('Number of Methylated Sites in Exons Across Samples')
plt.xlabel('Cancer Type')
plt.ylabel('# of methylated sites in exons per sample')
plt.tight_layout()
def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')

# autolabel(rects1)
# autolabel(rects2)


# In[5]:

a = 'download?fn=%2Frelease_23%2FProjects%2FPAEN-AU%2Fsimple_somatic_mutation.open.PAEN-AU.tsv'
cancer_type = a.split('.')[2]
print (cancer_type)


# In[ ]:



