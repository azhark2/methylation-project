
# coding: utf-8

# In[8]:

import numpy as np
import pickle as cPickle
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
actual_ratios_no_threshold = cPickle.load(open("/Users/khandekara2/Documents/methylationProject/03_results/actual_ratios_no_threshold.pickle", "rb"))


# In[9]:

plt.hist(actual_ratios_no_threshold, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Methylation Ratio')
plt.ylabel('Frequency of Overlaps')


# In[10]:

print len(actual_ratios_no_threshold)


# In[24]:

import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_MALY_cds.tsv', sep='\t')
samples = df['submitted_sample_id'].unique()
df['submitted_sample_id'] = df['submitted_sample_id'].str.strip()
df = df[df['submitted_sample_id'] == samples[0]]
df.info()
df['chromosome_strand'].unique()


# In[17]:

df.to_csv('CLLE_sample.csv')


# In[21]:

bed_files = []
prefix = 'CpG_hg19_chr'
suffix= '.bed'
for i in range(1, 23):
    bed_files.append(prefix + str(i) + suffix)
print bed_files


# In[22]:

x = 2.19e+08 
y = int(float(x))
print y


# In[23]:

type(x)


# In[ ]:



