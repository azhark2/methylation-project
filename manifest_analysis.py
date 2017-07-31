
# coding: utf-8

# In[12]:

import pandas as pd
manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
df = pd.read_csv(manifest_location, header=7)


# In[13]:

df.head


# In[16]:

df[['IlmnID', 'UCSC_RefGene_Name', 'CHR','MAPINFO', 'UCSC_RefGene_Group',  ]]


# In[8]:

print (len(genes))


# In[10]:

df = df.dropna(subset=['UCSC_RefGene_Name'])


# In[12]:

df.head(:, [1])


# In[13]:

df.iloc[5, [1, 2]]


# In[17]:

df = df.loc[df.iloc[7, 0]]
df.head()


# In[6]:

import pickle
dict2 = pickle.load( open( "master_dict.pickle", "rb" ) )


# In[20]:

print (dict.items(:5))


# In[22]:

print (dict2.items())


# In[10]:

for key in dict2:
    gene_tuple = dict2[key]
    gene = gene_tuple[0].split(';')
    dict2[key] = gene
    
print (dict2.items())


# In[17]:

df.head()


# In[25]:

df[(df.UCSC_RefGene_Name.string.contains('Exon'))]


# In[35]:

# df[(df['nationality'] == "USA")]
exons = df[(df['UCSC_RefGene_Group'].notnull()) & (df['UCSC_RefGene_Group'].str.contains('Exon'))]


# In[37]:

exons.head


# In[38]:

print (exons.shape[0])


# In[42]:

exon_ids = list(exons['IlmnID'])
print (exon_ids[0:5])
print (len(exon_ids))
exon_ids = set(exons['IlmnID'])
print (len(exon_ids))


# In[ ]:



