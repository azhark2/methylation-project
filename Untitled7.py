
# coding: utf-8

# In[1]:

#create violin plots 
import pandas as pd
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import seaborn as sns


# In[13]:

#create 4 separate dataframes for each prev_base and plot violinplots
clle = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/CLLE_prevBase.tsv', sep = '\t')
pbca = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/PBCA_prevBase.tsv', sep = '\t')
#filter out methylation values less than 0.3
clle = clle[clle['methylation_ratio'] > 0]
pbca = pbca[pbca['methylation_ratio'] > 0]


# In[11]:

pbca['prev_base'] = pbca['prev_base'].str.upper()
clle['prev_base'] = clle['prev_base'].str.upper()


# In[6]:

#break up dataframes into different motifs
A = clle[clle['prev_base'].str.upper() == "A"]
C = clle[clle['prev_base'].str.upper() == "C"]
T = clle[clle['prev_base'].str.upper() == "T"]
G = clle[clle['prev_base'].str.upper() == "G"]


# In[9]:

sns.violinplot([A.methylation_ratio], orient="v")
sns.violinplot([C.methylation_ratio], orient="v")
sns.violinplot([T.methylation_ratio], orient="v")
sns.violinplot([G.methylation_ratio], orient="v")


# In[12]:

ax = sns.violinplot(x=pbca.prev_base, y=pbca.methylation_ratio)


# In[ ]:



