
# coding: utf-8

# In[14]:

import csv
# import pickle
# list_of_coding_cpgs = set([]) #coordinates of coding cpg's
# bed_files = ['chr1_cds_cpg.bed', 'chr2_cds_cpg.bed', 'chr3_cds_cpg.bed', 'chr4_cds_cpg.bed', 'chr5_cds_cpg.bed', 'chr6_cds_cpg.bed', 'chr7_cds_cpg.bed',
#             'chr8_cds_cpg.bed', 'chr9_cds_cpg.bed', 'chr10_cds_cpg.bed', 'chr11_cds_cpg.bed', 'chr12_cds_cpg.bed', 'chr13_cds_cpg.bed', 'chr14_cds_cpg.bed',
#             'chr15_cds_cpg.bed', 'chr16_cds_cpg.bed', 'chr17_cds_cpg.bed', 'chr18_cds_cpg.bed', 'chr19_cds_cpg.bed', 'chr20_cds_cpg.bed',
#             'chr21_cds_cpg.bed', 'chr22_cds_cpg.bed']
# for file in bed_files:
#     with open(file, 'r') as f:
#         reader = csv.reader(f, delimiter='\t')
#         for row in reader:
#             list_of_coding_cpgs.add((row[0], row[1]))
# pickle.dump(list_of_coding_cpgs, open('list_of_coding_cpgs.pickle', 'wb'))

biseq_files = ['/Users/khandekara2/Documents/methylationProject/02_code/split_zz']
for file in biseq_files:
    with open(file, 'r') as f, open(file.split('.')[0] + '.bed', 'w') as csvout:
        reader = csv.reader(f, delimiter='\t')
        csvout = csv.writer(csvout, delimiter='\t')
        rownum = 0
        for row in reader: #take care of header
            location = ('chr' + str(row[6]), str(row[7]))
            if rownum == 0:
                csvout.writerow(row[3:14])
            elif location in list_of_coding_cpgs: #filter for only coding cpg's
                csvout.writerow(row[3:14])
            rownum += 1



        


# In[2]:

file = 'meth_seq_CLLE.tsv'
print(file.split('.')[0] + '.bed')


# In[ ]:

('chr1', '55066997')


# In[9]:

import pandas as pd
df = pd.read_csv('/Users/khandekara2/Documents/methylationProject/02_code/split_zz', sep='\t')


# In[10]:

df.head()


# In[11]:

df.info()


# In[16]:

x = 235537026
x = str(x)
print (x)


# In[17]:

type(x)


# In[19]:

s = 'download?fn=%2Frelease_23%2FProjects%2FPAEN-AU%2Fsimple_somatic_mutation.open.PAEN-AU.tsv'
split = s.split('.')
print (split[2] + '_mutation.bed')


# In[20]:

print (s.split('.')[2] + '_mutation.bed')


# In[ ]:



