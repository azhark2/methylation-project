
# coding: utf-8

# In[ ]:




# In[3]:

#create bed file of cds of cosmic cancer genes
import csv
import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/Census_allThu Feb  9 18_47_50 2017.csv', header=0)
df.head()
locations = df['Genome Location']
names = df['Gene Symbol']
with open('cancer_genes.bed', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    for name, location in zip(names, locations):
        loc = location.split(':')
        chrom = loc[0]
        start = loc[1].split('-')[0]
        stop =  loc[1].split('-')[1]
        if loc[1] != '-':
            writer.writerow(['chr' + chrom, start, stop, name])

        


# In[21]:

df.head()


# In[2]:

s = '10:26748570-26860863'
split = s.split(':')
print ('chr' + split[0])


# In[26]:

print (split[1].split('-')[0])
print (split[1].split('-')[1])


# In[ ]:

s1 = 'download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv'

