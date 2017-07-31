
# coding: utf-8

# In[114]:

import csv
import pandas as pd
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import seaborn as sns
pbca = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/MALY_prevBase_cancerGenes.tsv', sep = '\t')


# In[115]:

pbca.shape


# In[116]:

import numpy as np
thresholds = [x / 10.0 for x in range(1, 7)]
for threshold in thresholds:
    with open("PBCA_census_sample_motif_results" + str(threshold*10) +'.csv', 'w') as csvout:
        writer = csv.writer(csvout)
        writer.writerow(['Sample', 'Percentage Methylated', "ACG", 'TCG', 'CCG', 'GCG' ])
        for sample in pbca['id'].unique():
            sub_pbca = pbca[pbca['id'] == sample]
            total = sub_pbca.shape[0]
            sub_pbca = sub_pbca[sub_pbca['methylation_ratio'] > threshold]
            methylated = sub_pbca.shape[0]
            sub_pbca['prev_base'] = sub_pbca['prev_base'].str.upper()
            groupby = sub_pbca['methylation_ratio'].groupby(sub_pbca['prev_base'])
            df = groupby.agg([np.mean])
            df.reset_index(inplace=True)
            for base, ratio in zip(df['prev_base'], df['mean']):
                if base == 'A':
                    A_mean = ratio
                if base == 'T':
                    T_mean = ratio
                if base == 'C':
                    C_mean = ratio
                if base == 'G':
                    G_mean = ratio
            writer.writerow([sample, (methylated/total) * 100, A_mean, T_mean, C_mean, G_mean])
        


# In[84]:


print (pbca['id'].unique())


# In[117]:

import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots(1,1) 
x = np.arange(2,10,2)
ax.set_xticks(x)

# Set number of ticks for x-axis
ax.set_xticks(x)
x_ticks_labels = ['ACG','TCG','CCG','GCG']
# Set ticks labels for x-axis
ax.set_xticklabels(x_ticks_labels, fontsize=18)
ax.set_xlabel('Motif', fontsize=18)
ax.set_ylabel('Average Methylation Ratio', fontsize=18)



thresholds = [x / 10.0 for x in range(1, 7)]
for threshold in thresholds:
    df = pd.read_csv("PBCA_census_sample_motif_results" + str(threshold*10) +'.csv')
    A_mean = df['ACG'].mean()
    T_mean = df['TCG'].mean()
    C_mean = df['CCG'].mean()
    G_mean = df['GCG'].mean()
    means = [A_mean, T_mean, C_mean, G_mean]
    ax.plot(x, means, label=threshold)
plt.legend(loc='best')


# In[86]:

pbca_cancer = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/PBCA_prevBase_cancerGenes.bed', sep = '\t')
maly_cancer = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/MALY_prevBase_cancerGenes.bed', sep = '\t')


# In[ ]:

plt.hist(pbca_cancer['methylation_ratio'])


# In[ ]:

plt.hist(maly_cancer['methylation_ratio'])


# In[ ]:

grouped = pbca_cancer['methylation_ratio'].groupby(pbca_cancer['chromosome'], pbca['start'], pbca['stop']).agg(np.mean, np.std) 
# maly_cancer['methylation_ratio'.groupby = 


# In[105]:

df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/PBCA_prevBase_TP53cds.bed', sep = '\t')


# In[106]:

locations = []
for chrom, start, stop in zip(list(df['chromosome']), list(df['start']), list(df['stop'])):
    locations.append((chrom, start, stop))
df['location'] = locations


# In[107]:

grouped = df['methylation_ratio'].groupby(df['location'])


# In[108]:

groupby = grouped.agg([np.mean, np.std])
groupby.reset_index(inplace=True)


# In[102]:

groupby.head()


# In[109]:

groupby.shape


# In[111]:

y = groupby['mean']
y_err = groupby['std']
# start = 7565097
# x = []
# for loc in groupby['location']:
#     dist = loc[1] - start
#     x.append(dist)

x = [i for i in range(1, groupby.shape[0] + 1)]


# In[112]:

print (x)


# In[113]:

fig = plt.figure(figsize=(20, 2))
ax = fig.add_subplot(111)
ax.set_xticks(x)
ax.errorbar(x, y, yerr=y_err, fmt='o')
ax.set_xlabel('CpG #', fontsize=12)
ax.set_ylabel('Average Methylation Ratio', fontsize=12)


# In[97]:

print (len(x))


# In[ ]:



