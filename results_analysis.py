
# coding: utf-8

# In[3]:

import pickle as cPickle
actual_ratios = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/actual_ratios.pickle", "rb"))
cancer_overlaps = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/cancer_overlaps.pickle", "rb"))
WGBS_avgBeta = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/WGBS_avgBeta.pickle", "rb"))


# In[6]:

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
values = []
for location in cancer_overlaps['PBCA'][0]:
    values.append(WGBS_avgBeta['PBCA'][location])
plt.hist(values, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Methylation Ratio')
plt.ylabel('Frequency of Overlaps')
plt.title('PBCA')


# In[7]:

values = []
for location in cancer_overlaps['MALY'][0]:
    values.append(WGBS_avgBeta['MALY'][location])
plt.hist(values, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Methylation Ratio')
plt.ylabel('Frequency of Overlaps')
plt.title('MALY')


# In[8]:

plt.hist(actual_ratios, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Methylation Ratio')
plt.ylabel('Frequency of Overlaps')
# plt.title('MALY')


# In[9]:

import collections
counter=collections.Counter(cancer_overlaps['MALY'][0])
print(counter)


# In[10]:

counter=collections.Counter(cancer_overlaps['PBCA'][0])
print(counter)


# In[17]:

print len(actual_ratios)


# In[15]:

actual_ratios_no_threshold = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/actual_ratios_no_threshold.pickle", "rb"))


# In[16]:


plt.hist(actual_ratios_no_threshold, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Methylation Ratio')
plt.ylabel('Frequency of Overlaps')


# In[26]:

import pandas as pd 
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/CLLE_prevBase.tsv', sep = '\t')



# In[21]:

df.head()


# In[27]:

groupby = df['methylation_ratio'].groupby(df['prev_base'].str.upper())
df = groupby.agg([np.mean])
df.reset_index(inplace=True)


# In[28]:

df.head()


# In[33]:

df2 = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/PBCA_prevBase.tsv', sep = '\t')


# In[34]:

groupby2 = df2['methylation_ratio'].groupby(df2['prev_base'].str.upper())
df2 = groupby2.agg([np.mean])
df2.reset_index(inplace=True)


# In[35]:

df2


# In[38]:

import matplotlib.pyplot as plt

for i, group in groupby:
    plt.figure()
    group.hist()


# In[41]:

#create 4 separate dataframes for each prev_base and plot boxplot
clle = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/CLLE_prevBase.tsv', sep = '\t')
pbca = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/PBCA_prevBase.tsv', sep = '\t')
A = pbca[pbca['prev_base'].str.upper() == "A"]
C = pbca[pbca['prev_base'].str.upper() == "C"]
T = pbca[pbca['prev_base'].str.upper() == "T"]
G = pbca[pbca['prev_base'].str.upper() == "G"]


# In[42]:

print A.shape
print C.shape
print T.shape
print G.shape


# In[43]:

print pbca.shape


# In[46]:

pbca['prev_base'] = pbca['prev_base'].str.upper()


# In[47]:

#plot boxplots
pbca.boxplot(by='prev_base', column='methylation_ratio')


# In[48]:

A = clle[clle['prev_base'].str.upper() == "A"]
C = clle[clle['prev_base'].str.upper() == "C"]
T = clle[clle['prev_base'].str.upper() == "T"]
G = clle[clle['prev_base'].str.upper() == "G"]


# In[49]:

print A.shape
print C.shape
print T.shape
print G.shape


# In[53]:

clle['prev_base'] = clle['prev_base'].str.upper()


# In[54]:

clle.boxplot(by='prev_base', column='methylation_ratio')


# In[56]:


grouped = pbca['methylation_ratio'].groupby(pbca['id'])
grouped_df = grouped.agg([np.mean, np.std])
grouped_df.reset_index(inplace=True)


# In[57]:

grouped_df


# In[58]:

grouped_df.hist(column='mean')


# In[59]:

grouped_df.hist(column='std')


# In[4]:

import numpy as np
import pandas as pd
pbca = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/PBCA_prevBase.tsv', sep = '\t')
methylated = pbca[pbca['methylation_ratio'] > 0]
grouped = methylated['methylation_ratio'].groupby(methylated['id'])
grouped_df = grouped.agg([np.mean, np.std])
grouped_df.reset_index(inplace=True)


# In[6]:

get_ipython().magic('matplotlib inline')
grouped_df.hist(column='mean')


# In[7]:

grouped_df.hist(column='std')


# In[10]:

grouped_df.to_csv('sample_means.csv')


# In[ ]:



