
# coding: utf-8

# In[4]:

import os
import pandas as pd
os.chdir('/Users/khandekara2/Documents/methylationProject/01_data/hotspots')


# In[18]:

df = pd.read_csv('cancer_cds_cpgs_table.tsv', sep='\t')


# In[19]:

df.head()


# In[22]:

df.sort_values('Number of CpGs in CDS', inplace=True, ascending=False)


# In[28]:

df.head(20)


# In[31]:

get_ipython().magic('matplotlib inline')
df.hist(column='Number of CpGs in CDS', bins=20)


# In[29]:

tp53 = df[df['Gene'] == 'TP53']


# In[30]:

tp53.head()


# In[6]:

pbca = pd.read_csv('PBCA_expanded.tsv', sep='\t')
maly = pd.read_csv('MALY_expanded.tsv', sep='\t')
clle = pd.read_csv('CLLE_expanded.tsv', sep='\t')


# In[7]:

locations = []
for chrom, start, stop in zip(list(pbca['chromosome']), list(pbca['start']), list(pbca['stop'])):
    locations.append((chrom, start, stop))
pbca['location'] = locations

locations = []
for chrom, start, stop in zip(list(maly['chromosome']), list(maly['start']), list(maly['stop'])):
    locations.append((chrom, start, stop))
maly['location'] = locations


# In[8]:

import numpy as np
grouped = pbca['methylation_ratio'].groupby(pbca['location'])
pbca_grouped = grouped.agg([np.mean, np.std])
pbca_grouped.reset_index(inplace=True)

grouped = maly['methylation_ratio'].groupby(maly['location'])
maly_grouped = grouped.agg([np.mean, np.std])
maly_grouped.reset_index(inplace=True)



# In[11]:

maly_grouped.head()


# In[ ]:

pbca_avgBeta = {}
maly_avgBeta = {}
for location, mean, std in zip(pbca_grouped['location'], pbca_grouped['mean'], pbca_grouped['std']):
    pbca_avgBeta[str(location)] = (mean, std)

for location, mean, std in zip(maly_grouped['location'], maly_grouped['mean'], maly_grouped['std']):
    maly_avgBeta[str(location)] = (mean, std)


# In[15]:

print (maly_avgBeta.keys())


# In[16]:

df = pd.read_csv('all_hotspots.tsv', sep='\t')


# In[43]:

df.head()


# In[17]:

pbca = df[df['cancer_type'] == 'PBCA']
maly = df[df['cancer_type'] == 'MALY']


# In[23]:

maly.head()


# In[40]:

pbca_means = []
maly_means = []
for x in pbca['location']:
    location = x.replace("'", "")
    if location in pbca_avgBeta.keys():
        pbca_means.append(pbca_avgBeta[location])
    
for x in maly['location']:
    location = x.replace("'", "")
    if location in maly_avgBeta.keys():
        maly_means.append(maly_avgBeta[location])


# In[26]:

print (maly_avgBeta[('chr11', '64572568', '64572569)])


# In[41]:

print (maly_means)


# In[31]:

print (('chr11', 64572568, 64572569) in maly_avgBeta.keys())


# In[34]:

for location in pbca['location']:
    print (type(location[0]))
    print (type(location[1]))
    print (type(location[2]))
    print (type(location))
    break


# In[ ]:



