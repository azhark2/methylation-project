
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
import pandas as pd
import os
os.chdir('/Users/khandekara2/Documents/methylationProject/01_data/hotspots')
df = pd.read_csv('all_hotspots.tsv', sep='\t')


# In[2]:

df.info()


# In[3]:

df.head()


# In[5]:

df2 = pd.read_csv('gene_id_names.tsv', sep='\t')


# In[7]:

df2.head()


# In[14]:

#create dictionary that maps each Ensembl id to a tuple with (gene name, gene description)
idToGene = {}
for ids, name, description in zip(df2['Gene stable ID'], df2['Gene name'], df2['Gene description']):
    idToGene[ids] = (name, description)


# In[16]:

genes = []
descriptions = []
missing = 0
for ids in df['gene_id']:
    if ids in idToGene:
        gene = idToGene[ids][0]
        description = idToGene[ids][1]
        genes.append(gene)
        descriptions.append(description)
    else:
        missing += 1
df['gene'] = genes
df['description'] = descriptions
print (missing)


# In[15]:

print (idToGene)


# In[17]:

df.head()


# In[18]:

df['gene'].value_counts()


# In[19]:

locations = []
for chrom, start, stop in zip(df['chromosome'], df['start'], df['stop']):
    locations.append((chrom, start, stop))
df['location'] = locations
df.head()


# In[20]:

df['location'].value_counts()


# In[21]:

smarca = df[df['gene'] == 'SMARCA4']


# In[23]:

smarca.head(7)


# In[24]:

df.head()


# In[25]:

df.to_csv('all_hotspots.tsv', sep='\t', index=False)


# In[26]:

smarca.to_csv('smarca_hotspots.tsv', sep='\t', index=False)


# In[27]:

tp53 = df[df['gene'] == 'TP53']


# In[29]:

tp53.to_csv('tp53_hotspots.tsv', sep='\t', index=False)


# In[30]:

tp53.head()


# In[31]:

df.info()


# In[33]:

cancer_types = []
for cancer in df['id']:
    if cancer.startswith('ICGC'):
        cancer_types.append('PBCA')
    elif cancer.startswith('tumor'):
        cancer_types.append('MALY')
    else:
        cancer_types.append('?')
print (cancer_types)


# In[34]:

df['cancer_type'] = cancer_types


# In[37]:

df.to_csv('all_hotspots.tsv', sep='\t', index=False)
df.to_csv('all_hotspots.bed', sep='\t', index=False, header=False)


# In[ ]:

pbca = df[df['cancer_type'] == 'PBCA']
maly = df[df['cancer_type'] == 'MALY']

