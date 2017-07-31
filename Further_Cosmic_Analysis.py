
# coding: utf-8

# In[1]:

import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/CosmicCompleteDifferentialMethylation.csv', header = None)


# In[2]:

df.head()


# In[5]:

df.columns = ['study_id', 'id_sample', 'sample_name', 'id_tumour', 'primary_site', 'site_subtype1', 'site_subtype2', 'site_subtype3', 'primary_histology',
              'histology_subtype1', 'histology_subtype2', 'histology_subtype3', 'fragment_id', 'genome_version', 'chr', 'position', 'gene_name',
              'methylation', 'avg_beta_value_normal', 'beta_value', 'p_value',]


# In[6]:

non_promoter = df[df.gene_name.str.contains('Promoter') == False]
print len(non_promoter)


# In[7]:

print df.shape


# In[8]:

promoter = df[df.gene_name.str.contains('Promoter') == True]


# In[9]:

print len(promoter)


# In[10]:

print len(promoter) + len(non_promoter)


# In[11]:


print promoter.head()





# In[12]:

intragenic = df[(df.gene_name.str.contains('Unclassified') == False) & (df.gene_name.str.contains('Promoter') == False)
                & (df.gene_name.str.contains('Associated') == False)]


# In[13]:

intragenic.head()


# In[14]:

df.describe()


# In[15]:

df2 = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/Census_allThu Feb  9 18_47_50 2017.csv', )


# In[16]:

df2 = df2.fillna('')


# df2.shape
# 
# 

# In[17]:

print df2.shape


# In[ ]:




# In[ ]:




# In[18]:

genes = df2.iloc[0]


# In[19]:

print len(genes)


# In[20]:

genes = df2.iloc[:, 0]


# In[21]:

print len(genes)


# In[22]:

gene_list = list(genes)


# In[23]:

print gene_list


# In[24]:

gene_type = df2.iloc[:, 12]


# In[25]:

gene_type = list(gene_type)


# In[26]:

print gene_type


# In[27]:

result = zip(gene_list, gene_type)


# In[28]:

gene_dict = {}
for g in result:
    gene_dict[g[0]] = g[1]
print gene_dict
    


# In[29]:

intragenics = list(intragenic)
print len(intragenics)


# In[ ]:




# In[30]:

intragenic = df[(df.gene_name.str.contains('Unclassified') == False) & (df.gene_name.str.contains('Promoter') == False)
                & (df.gene_name.str.contains('Associated') == False)]
print len(intragenics)


# In[31]:

print intragenics


# In[32]:

print len(intragenic)


# In[33]:

intragenics = list(intragenic)
count = 0
for gene in intragenics:
    if gene in gene_dict:
        count +=1
print count


# In[34]:

print intragenics


# In[35]:

print intragenic


# In[36]:

intragenic = df[(df.gene_name.str.contains('Unclassified') == False) & (df.gene_name.str.contains('Promoter') == False)
                & (df.gene_name.str.contains('Associated') == False)]
print len(intragenic)


# In[37]:

print list(intragenic)


# In[38]:

intragenics = df.iloc[:, 16]
print len(intragenics)


# In[65]:

intragenic = set(intragenics)
count = 0
nan_count = 0
ts = 0
onco = 0
both = 0
for gene in intragenic:
    if gene in gene_dict and gene_dict[gene]!= '':
        count +=1
        if gene_dict[gene] == 'TSG':
            ts += 1
        elif gene_dict[gene] == 'oncogene':
            onco += 1
        else:
            both += 1     
    elif gene in gene_dict and gene_dict[gene] == '':
        nan_count += 1
print count
print nan_count
print ts
print onco
print both
print ts + onco + both


# In[66]:

print 'There are %d differentially methylated sites that fall in the coding regions of  tumor suppressors or oncogenes' %count
print "%d of these are tumor suppressors" % ts
print "%d of these are oncogenes" % onco
print "%d of these are annotated as both oncogenes and tumor suppressors?" % both


# In[41]:

nan = gene_dict['TOP1']
print nan


# In[ ]:




# In[42]:

positions = list(df.iloc[:, 15])
print len(positions)


# In[43]:

breast_positions = df[df.primary_site.str.contains('breast') == True]
bladder_positions = df[df.primary_site.str.contains('urinary') == True]


# In[44]:

print len(breast_positions)
print len(bladder_positions)


# In[45]:

breast = list(breast_positions.iloc[:, 15])
bladder = list(bladder_positions.iloc[:, 15])
print len(breast_positions)
print len(bladder_positions)


# In[46]:

print breast
print bladder


# In[47]:

print gene_dict['TP53']


# In[48]:

p53 = df[df.position >= 0 and df.position < 100000]


# In[ ]:

p53 = df[df.position > 0]


# In[ ]:

print len(p53)


# In[ ]:

chr17 = df[df.chr == 17]


# In[ ]:

print len(chr17)


# In[49]:

p53 = chr17[(chr17.position >= 7668402) & (chr17.position < 7687538)]


# In[50]:

print len(p53)



# In[51]:

if 'TP53' in list(intragenic):
    print 'Yes'
    


# In[ ]:




# In[67]:

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
labels = 'Tumor Suppressors', 'Oncogenes', 'TSG/Oncogene'
sizes = [ts, onco, both]
explode = (0, 0, 0)
fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.show()


# In[55]:

df.gene_name.values


# In[56]:

df.gene_name.values.unique()


# In[57]:

df[gene_name].values.unique()


# In[60]:

df[df.gene_name].nunique(dropna=True)


# In[64]:

gene = df[df.gene_name.dropna()]


# In[ ]:

df3 = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/raw_data/ICGC/meth_array.PAEN-AU.tsv', sep='\t', header = None)


# In[ ]:



