
# coding: utf-8

# In[10]:

#write out ICGC somatic mutation files as bed files
import pandas as pd

all_ssms = ['download?fn=%2Frelease_23%2FProjects%2FPAEN-AU%2Fsimple_somatic_mutation.open.PAEN-AU.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBLCA-US%2Fsimple_somatic_mutation.open.BLCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBRCA-US%2Fsimple_somatic_mutation.open.BRCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCESC-US%2Fsimple_somatic_mutation.open.CESC-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCOAD-US%2Fsimple_somatic_mutation.open.COAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FGBM-US%2Fsimple_somatic_mutation.open.GBM-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FKIRC-US%2Fsimple_somatic_mutation.open.KIRC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FKIRP-US%2Fsimple_somatic_mutation.open.KIRP-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FLGG-US%2Fsimple_somatic_mutation.open.LGG-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FLIHC-US%2Fsimple_somatic_mutation.open.LIHC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FOV-US%2Fsimple_somatic_mutation.open.OV-US.tsv',
         'download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv', 'download?fn=%2Frelease_23%2FProjects%2FREAD-US%2Fsimple_somatic_mutation.open.READ-US.tsv'
        'download?fn=%2Frelease_23%2FProjects%2FSKCM-US%2Fsimple_somatic_mutation.open.SKCM-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FSTAD-US%2Fsimple_somatic_mutation.open.STAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FTHCA-US%2Fsimple_somatic_mutation.open.THCA-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FUCEC-US%2Fsimple_somatic_mutation.open.UCEC-US.tsv']

ssms = ['/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', '/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Fcurrent%2FProjects%2FMALY-DE%2Fsimple_somatic_mutation.open.MALY-DE.tsv', '/Users/khandekara2/Documents/methylationProject/01_data/download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv']

for file in ssms:
    df = pd.read_csv(file, sep='\t')
    df1 = df[(df.mutated_from_allele == 'C') & (df.mutated_to_allele == 'T')]
    df2 = df[(df.mutated_from_allele == 'G') & (df.mutated_to_allele == 'A')]
    df = pd.concat([df1, df2])
    cols = list(df.columns) 
#     df = df[['mean', '0', '1', '2', '3']]


# In[11]:

print (cols)


# In[12]:

df = df[['chromosome', 'chromosome_start', 'chromosome_end', 'submitted_sample_id', 'mutated_from_allele', 'mutated_to_allele', 'gene_affected']]


# In[13]:

df.to_csv(file.split('.')[2] + '_mutation.bed')


# In[14]:

f = file.split('.')


# In[24]:

name = file.split('.')[2] + '_mutation.bed'


# In[25]:

print (name)


# In[ ]:



