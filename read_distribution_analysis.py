
# coding: utf-8

# # Analysis of the read depth distributions of the bisulfite sequencing data from ICGC. The cancers analyzed are PBCA, BOCA, CLLE, and MALY. Coverage represents the number of times a sequenced DNA fragment (i.e., a read) maps to a genomic target. The deeper the coverage of a target region (i.e., the more times the region is sequenced), the greater the reliability and sensitivity of the sequencing assay. Typically, the minimum depth of coverage required for genomic resequencing of diploid organisms, such as human, mouse, or rat, is 20–30X. 

# In[16]:

import pandas as pd
pbca = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_PBCA_cds.tsv', sep='\t')
boca = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_BOCA_cds.tsv', sep='\t')
clle = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_CLLE_cds.tsv', sep='\t')
maly = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_MALY_cds.tsv', sep='\t')


# In[17]:

#make new column that contains total read count(methylated + unmethylated)
pbca_totalReads = []
clle_totalReads = []
boca_totalReads = []
maly_totalReads = []

for m, u in zip(pbca['methylated_reads'], pbca['unmethylated_reads']):
    pbca_totalReads.append(int(m + u))

for m, u in zip(clle['methylated_reads'], clle['unmethylated_reads']):
    clle_totalReads.append(int(m + u))
    
for m, u in zip(boca['methylated_reads'], boca['unmethylated_reads']):
    boca_totalReads.append(int(m + u))
    
for m, u in zip(maly['methylated_reads'], maly['unmethylated_reads']):
    maly_totalReads.append(int(m + u))


# In[18]:

pbca['total_reads'] = pbca_totalReads
clle['total_reads'] = clle_totalReads
boca['total_reads'] = boca_totalReads
maly['total_reads'] = maly_totalReads


# In[19]:

get_ipython().magic('matplotlib inline')
pbca.hist(column ='total_reads', bins = 100, range=[0, 100])


# In[20]:

clle.hist(column = 'total_reads', bins = 100, range=[0, 100])


# In[21]:

maly.hist(column = 'total_reads', bins = 100, range=[0, 100])


# In[22]:

boca.hist(column = 'total_reads', bins = 100, range=[0, 100])


# In[23]:

pbca_underTen = 0
clle_underTen = 0
maly_underTen = 0
for pb, cl, ma in zip(pbca['total_reads'], clle['total_reads'], maly['total_reads']):
    if pb < 10:
        pbca_underTen += 1
    if cl < 10:
        clle_underTen += 1
    if ma < 10:
        maly_underTen += 1
print (pbca_underTen)
print (clle_underTen)


# In[24]:

#calculate mean mapped read depth(total # of reads / total of bases) for each cancer
pbca.info()


# In[ ]:



