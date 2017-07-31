
# coding: utf-8

# In[9]:

import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/CosmicCompleteDifferentialMethylation.csv', header = None)


# In[10]:

df.head()


# In[46]:

df.columns = ['study_id', 'id_sample', 'sample_name', 'id_tumour', 'primary_site', 'site_subtype1', 'site_subtype2', 'site_subtype3', 'primary_histology',
              'histology_subtype1', 'histology_subtype2', 'histology_subtype3', 'fragment_id', 'genome_version', 'chr', 'position', 'gene_name',
              'methylation', 'avg_beta_value_normal', 'beta_value', 'p_value',]


# In[12]:

promoter = df[df.gene_name.str.contains('Promoter') == True]
print len(promoter)


# In[13]:

gene = df[df.gene_name.dropna()]


# In[38]:

unclassified = df[df.gene_name.str.contains('Unclassified') == True]
print len(unclassified)


# In[18]:

gene_associated = df[df.gene_name.str.contains('Gene_Associated') == True]
print len(gene_associated)


# In[24]:

nulls = df.isnull().values.sum()
print nulls


# In[31]:

coding = df[(df.gene_name.str.contains('Unclassified') == False) & (df.gene_name.str.contains('Promoter') == False)
                & (df.gene_name.str.contains('Associated') == False) & (df.gene_name.str.contains('NaN') == False)]
print len(coding)


# In[40]:

print len(coding) + nulls + len(unclassified) + len(promoter) + len(gene_associated)


# In[27]:

df.shape


# In[32]:

sum(pd.isnull(df['gene_name']))


# In[44]:

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
labels = 'Unclassified', 'Coding Region', 'Missing', 'Promoter/Gene Associated'
sizes = [len(unclassified), len(coding), nulls, len(promoter) + len(gene_associated)]
explode = (0, 0, 0, 0)
fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.show()


# In[49]:

print len(df['id_sample'].unique())


# In[50]:

print len(df['id_tumour'].unique())


# In[56]:

positions = list(df['position'])


# In[53]:

print positions

    


# In[57]:

print len(positions)
positions_set = set(positions)
print len(positions_set)


# In[69]:

# first ensure there are no two identical CpG sites on different chromosomes
import collections
counter=collections.Counter(positions)
repeat_positions = []
for key in counter:
    if counter[key] >= 465:
        repeat_positions.append(key)
print len(repeat_positions)



# In[67]:

print max(counter.values())


# In[70]:

for key in counter.keys():
    if counter[key] == 464:
        print key

    


# In[1]:

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
list = [2, 2, 1, 1, 1, 1, 3, 3, 3]
plt.hist(list, bins = 3)
plt.show()


# In[6]:

filename = 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PAEN-AU/meth_array.PAEN-AU.tsv.gz'
f = 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PAEN-AU/meth_array.PAEN-AU.tsv.gz'
print("finished reading and pickling " + f[len(f) - 14:-6])


# In[14]:


numMethylated = 240
print ('Number of methylated sites for ' + f[len(f) - 14:-7] + ': %d' %(numMethylated))   
print ('Number of total methylated sites across all cancer types: %d ' %(numMethylated))
print ('Number of total methylated sites across all cancer types that fall into the coding regions of genes: %d' %(numMethylated))


# In[ ]:



