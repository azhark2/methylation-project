
# coding: utf-8

# In[1]:

import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_CLLE_cds.tsv', sep='\t')


# In[2]:

df.head()


# In[3]:

df.shape


# In[4]:

df['submitted_sample_id'].value_counts()


# In[5]:

df['icgc_sample_id'].value_counts()


# In[6]:

df['analysis_id'].value_counts()


# In[7]:

df.info()


# In[8]:

methylated = df[df['methylation_ratio'] > 0.8]


# In[9]:

methylated.shape


# In[11]:

import csv
cancer_cpgs = set([])
file = '/Users/khandekara2/Documents/methylationProject/01_data/cancer_cds_cpg.bed'
with open(file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        cancer_cpgs.add((row[0], row[1]))
        
 #dump set of all CpG sites that are in cds region(based on GRCH38) in pickle file


# In[12]:

print (cancer_cpgs.pop())


# In[14]:

count = 0

for chrom, coord in zip(list(methylated['chromosome']), list(methylated['chromosome_start'])):
    location = ('chr' + str(chrom), str(coord))
    if location in cancer_cpgs:
        count += 1
print (count)


# In[15]:

count = 0

for chrom, coord in zip(list(df['chromosome']), list(df['chromosome_start'])):
    location = ('chr' + str(chrom), str(coord))
    if location in cancer_cpgs:
        count += 1
print (count)


# In[18]:

ssm = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', sep='\t')


# In[19]:

ssm = ssm[
        (ssm.reference_genome_allele == 'C') & (ssm.mutated_to_allele == 'T') | (ssm.reference_genome_allele == 'G') & (
        ssm.mutated_to_allele == 'A')]


# In[20]:

ssm = ssm[ssm.consequence_type == 'missense_variant']


# In[21]:

ssm.shape


# In[22]:

df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', sep='\t')


# In[23]:

df['submitted_sample_id'].value_counts()


# In[25]:

df['icgc_sample_id'].unique()


# In[26]:

print(len(df['icgc_sample_id'].unique()))


# In[30]:

count = 0

for chrom, coord in zip(list(ssm['chromosome']), list(ssm['chromosome_start'])):
    location = ('chr' + str(chrom), str(coord))
    if location in cancer_cpgs:
        count += 1
print (count)


# In[90]:

#Goal: Figue out if there are overlapping samples betwen bisulfite sequencing data and somatic mutation data

import pandas as pd

    
biseq = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_CLLE_cds.tsv', sep='\t')
biseq_samples = set(biseq['submitted_sample_id'].unique())
mutation = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Fcurrent%2FProjects%2FMALY-DE%2Fsimple_somatic_mutation.open.MALY-DE.tsv', sep='\t')
mutation_samples = set(mutation['submitted_sample_id'].unique())
common_samples = biseq_samples & mutation_samples
print (len(biseq_samples))
print (len(mutation_samples))
print (len(common_samples))


# In[ ]:




# In[6]:

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
s = (2, 3, 1)
#     2,Ab
#     3,aB
#     1,AB
v = venn2(subsets=s, set_labels=('WGBS Samples', 'Somatic Mutation Samples'))
v.get_label_by_id('10').set_text('26')
v.get_label_by_id('01').set_text('100')
v.get_label_by_id('11').set_text('22')
v.get_patch_by_id('10').set_color('c')
v.get_patch_by_id('01').set_color('#993333')
v.get_patch_by_id('11').set_color('blue')
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)
c = venn2_circles(subsets=s, linestyle='solid')
c[0].set_ls('dashed')  # Line style
c[0].set_lw(2.0)       # Line width
plt.show()


# In[46]:

import csv
cancer_cpgs = set([])
file = '/Users/khandekara2/Documents/methylationProject/01_data/cancer_cds_cpg.bed'
with open(file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        cancer_cpgs.add((row[0], (row[1])))
    print (type(row[0]))
    print (type(row[1]))


# In[48]:

#Goal: to do a sample by sample comparison of mutated and methylated sites CpG sites and write them out to csv file

        
with open("sampletosample_results.csv", 'w') as csvout:
    writer = csv.writer(csvout)
    #header
    for sample in list(common_samples):
        #subset dataframes for specific samples
        sub_biseq = biseq[biseq['submitted_sample_id'].str.contains(sample)]
        sub_mutation = mutation[mutation['submitted_sample_id'].str.contains(sample)]
        methylated = sub_biseq[sub_biseq['methylation_ratio'] > 0.7]
        methylated_cds = set([]) #all methylated sites in cds
        methylated_cancer_cds_cpg = set([]) #all methylated sites in cancer cds
        mutation_cds = set([])
        mutation_cancer_cds = set([]) # all mutated sites in cancer cds cpgs
        
        #find methylated sites in 
        for chrom, coord in zip(list(methylated['chromosome']), list(methylated['chromosome_start'])):
            location = ('chr' + str(chrom), str(coord))
            methylated_cds.add(location)
            if location in cancer_cpgs:
                methylated_cancer_cds_cpg.add(location)
        
        
        sub_mutation = sub_mutation[
        (sub_mutation.reference_genome_allele == 'C') & (sub_mutation.mutated_to_allele == 'T') | (sub_mutation.reference_genome_allele == 'G') & (
        sub_mutation.mutated_to_allele == 'A')]
        sub_mutation= sub_mutation[sub_mutation.consequence_type == 'missense_variant']
        
        asm = sub_biseq[sub_biseq['methylation_ratio'] > 0.4] # allele specific methylation
        asm_cds = set([]) # all partially methylated sites
        asm_cancer_cds = set([])
        
        #find partially methylated sites in cds and cancer cds
        for chrom, coord in zip(list(asm['chromosome']), list(asm['chromosome_start'])):
            location = ('chr' + str(chrom), str(coord))
            asm_cds.add(location)
            if location in cancer_cpgs:
                asm_cancer_cds.add(location)
        
        #find mutated sites in cds and cancer cds     
        for chrom, coord in zip(list(sub_mutation['chromosome']), list(sub_mutation['chromosome_start'])):
            location = ('chr' + str(chrom), str(coord))
            mutation_cds.add(location)
            if location in cancer_cpgs:
                mutation_cancer_cds.add(location) 
        
        #lists of overlaps 
        cds_overlap = asm_cds and mutation_cds #mutated and methylated sites in cds
        cancer_cds_overlap = asm_cancer_cds and mutation_cancer_cds ##mutated and methylated sites in cancer cds 
        
        #write results to csv file:
        #num sites covered, num cds sites methylated, num cancer cds sites methylated, num non cancer cds sites methylated,
        #num cds cpg sites mutated, num cancer cds cpg sites mutated, num non cancer cds cpg sites mutated,     
        writer.writerow([sample, sub_biseq.shape[0], len(methylated_cds), len(methylated_cancer_cds_cpg), len(methylated_cds) - 
                        len(methylated_cancer_cds_cpg), sub_mutation.shape[0], len(mutation_cancer_cds),
                        sub_mutation.shape[0] - len(mutation_cancer_cds), len(cds_overlap), len(cancer_cds_overlap),
                        len(cds_overlap) - len(cancer_cds_overlap)]) 
            
            

        
        
        
   


# In[91]:

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
biseq.hist(column='methylation_ratio', bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))


# In[57]:

import numpy as np
avgBeta = {}
locations = []
for chrom, coord in zip(list(biseq['chromosome']), list(biseq['chromosome_start'])):
    locations.append(('chr' + str(chrom), str(coord)))
biseq['location'] = locations
# new = biseq.filter(['location', 'methylation_ratio'], axis=1)
groupby = biseq['methylation_ratio'].groupby(biseq['location'])
df = groupby.agg([np.mean])
df.head()


# In[54]:

df.reset_index(inplace=True) 


# In[55]:

df.head()


# In[56]:

cancer_type = wgbs_file.split('_')[2]
for location, mean in zip(df['location'], df['mean']):
    avgBeta[cancer_type][location] = avg


# In[58]:

df.shape


# In[87]:

biseq = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_MALY_cds.tsv', sep='\t')
biseq_samples = set(biseq['submitted_sample_id'].unique())
mutation = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', sep='\t')
mutation_samples = set(mutation['submitted_sample_id'].unique())
common_samples = biseq_samples & mutation_samples
print (len(biseq_samples))
print (len(mutation_samples))
print (len(common_samples))


# In[61]:

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
s = (2, 3, 1)
#     2,Ab
#     3,aB
#     1,AB
v = venn2(subsets=s, set_labels=('WGBS Samples', 'Somatic Mutation Samples'))
v.get_label_by_id('10').set_text('5')
v.get_label_by_id('01').set_text('224')
v.get_label_by_id('11').set_text('2')
v.get_patch_by_id('10').set_color('c')
v.get_patch_by_id('01').set_color('#993333')
v.get_patch_by_id('11').set_color('blue')
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)
c = venn2_circles(subsets=s, linestyle='solid')
c[0].set_ls('dashed')  # Line style
c[0].set_lw(2.0)       # Line width
plt.show()


# In[79]:

biseq = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/meth_seq_PBCA_cds.tsv', sep='\t')
biseq_samples = set(biseq['submitted_sample_id'].unique())
mutation = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv', sep='\t')
mutation_samples = set(mutation['submitted_sample_id'].unique())
common_samples = biseq_samples & mutation_samples
print (len(biseq_samples))
print (len(mutation_samples))
print (len(common_samples))


# In[63]:

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
s = (2, 3, 1)
#     2,Ab
#     3,aB
#     1,AB
v = venn2(subsets=s, set_labels=('WGBS Samples', 'Somatic Mutation Samples'))
v.get_label_by_id('10').set_text('42')
v.get_label_by_id('01').set_text('380')
v.get_label_by_id('11').set_text('34')
v.get_patch_by_id('10').set_color('c')
v.get_patch_by_id('01').set_color('#993333')
v.get_patch_by_id('11').set_color('blue')
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)
c = venn2_circles(subsets=s, linestyle='solid')
c[0].set_ls('dashed')  # Line style
c[0].set_lw(2.0)       # Line width
plt.show()


# In[64]:

mutation.head()


# In[80]:

mut = mutation.iloc[:, 6:17]


# In[81]:

mut.head()


# In[ ]:



