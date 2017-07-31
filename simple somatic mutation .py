
# coding: utf-8

# In[8]:

import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/02_code/download?fn=%2Fcurrent%2FProjects%2FCOAD-US%2Fsimple_somatic_mutation.open.COAD-US.tsv', sep='\t')


# In[11]:

df = df[(df.reference_genome_allele == 'C') & (df.mutated_to_allele == 'T') | (df.reference_genome_allele == 'G') & (df.mutated_to_allele == 'A')]


# In[12]:

df.head()


# In[107]:

df.shape[0]


# In[15]:

df.mutated_to_allele


# In[16]:

df.info()


# In[19]:

df = df[df.consequence_type == 'missense_variant']


# In[20]:

df.info()


# In[30]:

chromosomes = df['chromosome']
positions = df['chromosome_start']
missense_mutations = []
for chromosome, position in zip(chromosomes, positions):
    missense_mutations.append((chromosome, position))

print (len(missense_mutations))
    


# In[31]:

df2 = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/raw_data/ICGC/HumanMethylation450_15017482_v1-2.csv', sep='/t')


# In[32]:

df.head()


# In[33]:

df2.head()


# In[21]:

manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
df2 = pd.read_csv(manifest_location, header=7)


# In[22]:

#create list of tuples of all assayed cpg site locations
cpg_sites = []
for c, p in zip(list(df2['CHR']), list(df2['MAPINFO'])):
    cpg_sites.append((c,p))
    


# In[39]:

overlaps = list(set(missense_mutations) & set(methylated_sites))


# In[41]:

print (len(overlaps))


# In[43]:

print (overlaps)


# In[48]:

print (len(df.icgc_donor_id.unique()))


# In[49]:

print (len(df.icgc_sample_id.unique()))


# In[ ]:

df3 = pd.read_csv('/Users/khandekara2/Documents/methylationProject/02_code/download?fn=%2Fcurrent%2FProjects%2FCOAD-US%2Fmeth_array.COAD-US.tsv', sep )


# In[3]:

import pandas as pd
tp = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/raw_data/ICGC/meth_array.PAEN-AU.tsv', sep='\t', iterator=True, chunksize=1000)
array = pd.concat(tp, ignore_index=True)
array.head()


# ssm = pd.read_csv('simple_somatic_mutation.open.PAEN-AU.tsv.gz', delimiter='/t')

# In[8]:

ssm = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/raw_data/ICGC/simple_somatic_mutation.open.PAEN-AU.tsv', sep='\t')


# In[11]:

ssm = ssm[(ssm.reference_genome_allele == 'C') & (ssm.mutated_to_allele == 'T') | (ssm.reference_genome_allele == 'G') & (ssm.mutated_to_allele == 'A')]
ssm.head()


# In[19]:

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
ssm.consequence_type.value_counts().plot(kind='barh')
plt.xscale('log')
plt.title('Number of appearances in dataset')
plt.xlabel('Frequency')


# In[20]:

ssm.shape


# In[23]:

#create list of tuples of all locations where C->T/G->A transitions occurred
chromosomes = ssm['chromosome']
positions = ssm['chromosome_start']
mutations = []
for chromosome, position in zip(chromosomes, positions):
    mutations.append((chromosome, position))



# In[27]:

overlaps = list(set(mutations) & set(cpg_sites))
print (len(overlaps))
print (overlaps)


# In[28]:

manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
manifest = pd.read_csv(manifest_location, header=7)


# In[29]:

manifest.UCSC_RefGene_Group.unique()


# In[47]:

master_dict = {}
gene_regions = manifest.UCSC_RefGene_Name
for c, p in zip(list(manifest['CHR'].dropna()), list(manifest['MAPINFO'].dropna())):
    master_dict[c, int(p)] = [0, 0]


# In[48]:

print (master_dict.items())


# In[42]:

print (len(list(manifest['CHR'].dropna())))


# In[46]:

probeToLocation = {}
for id, chr, coord in zip(list(manifest.IlmnID.dropna()), list(manifest['CHR'].dropna()), list(manifest['MAPINFO'].dropna())):
    probeToLocation[id] = (chr, int(coord))
print (probeToLocation.items())
# for chr, coord, beta in zip(array.)


# In[50]:

array.memory_usage(index=True).sum()


# In[ ]:


for probe, beta in zip(array.probe_id, array.methylation_value):
    master_dict[probeToLocation[probe]][0] += beta
    master_dict[probeToLocation[probe]][1] += 1
    


# In[51]:

array.submitted_sample_id.nunique()


# In[52]:

ssm.submitted_sample_id.nunique()


# In[53]:

print (ssm.submitted_sample_id.unique())


# In[55]:

print (array.submitted_sample_id.unique())


# In[65]:

array_samples = list(array.submitted_sample_id.unique())
ssm_samples = list(ssm.submitted_sample_id.unique())


# In[66]:

sample_overlaps = list(set(array_samples) & set(ssm_samples))
print(len(sample_overlaps))
print (sample_overlaps)


# In[60]:

list1 = [1, 2, 3]
list2 = [1, 2, 3, 4]
list(set(list1) & set(list2))


# In[63]:

type(ssm_samples[1])


# In[68]:

for sample in sample_overlaps:
    sample_array = array[array.submitted_sample_id == sample]
    sample_ssms = ssm[ssm.submitted_sample_id == sample]
    chromosomes = ssm['chromosome']
    positions = ssm['chromosome_start']
    types = ssm['consequence_type']
    chrs = array[]
    coords = array[]
    mutations = []
    cpg_sites = []
    for chromosome, position in zip(chromosomes, positions):
        mutations.append((chromosome, position))
    for chromosome, position, in zip (chrs, coords):
        cpg_sites.append((chromosome, position))
    #now look for overlaps
    overlaps = list(set(mutations) & set(cpg_sites))


# In[72]:

sample_array.submitted_sample_id.nunique()


# In[71]:

sample_ssms.submitted_sample_id.nunique()


# In[87]:

groupby = array['methylation_value'].groupby(array['probe_id'])


# In[81]:

type(groupby.mean())


# In[99]:

import numpy as np
df = groupby.agg([np.mean, np.std])
df.reset_index(inplace=True)


# In[100]:

df.info()


# In[105]:

cancertypeToAvgBeta = {}
cancertypeToAvgBeta['cancer_type'] = {}
for cg, avg, std in zip(list(df['probe_id']), list(df['mean']), list(df['std'])):
    cancertypeToAvgBeta['cancer_type'][cg] = (avg, std)


# In[106]:

cancertypeToAvgBeta['cancer_type'].items()


# In[ ]:



