
# coding: utf-8

# In[1]:

import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/raw_data/ICGC/HumanMethylation450_15017482_v1-2.csv', header = 7)


# In[2]:

df.head()


# In[3]:

df.info()
df.shape


# In[4]:

df.DMR.unique()


# In[6]:

print (len(df.IlmnID.unique()))


# In[7]:

df = df.dropna(subset = ['UCSC_RefGene_Name'])
df.shape


# In[9]:

genes = df.UCSC_RefGene_Name
print (len(genes))


# In[10]:

cg_ids = df.IlmnID 
print (len(cg_ids))
    


# In[11]:

id_dict = {}
for cg_id, gene in zip(cg_ids, genes):
    id_dict[cg_id] = gene

    


# In[17]:

count = 0
for key in id_dict:
    if count < 50:
        print (key)
        print (id_dict[key])
    count += 1
print (len(id_dict))


# In[ ]:



