
# coding: utf-8

# In[2]:

import pandas as pd
manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
df = pd.read_csv(manifest_location, header=7)
df = df.dropna(subset=['UCSC_RefGene_Name'])
exons = df[(df['UCSC_RefGene_Group'].notnull()) & (df['UCSC_RefGene_Group'].str.contains('Exon'))]
exon_ids = set(exons['IlmnID'])


# In[11]:

filenames = ['download?fn=%2Fcurrent%2FProjects%2FBLCA-US%2Fmeth_array.BLCA-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FBRCA-US%2Fmeth_array.BRCA-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FCESC-US%2Fmeth_array.CESC-US.tsv',
          'download?fn=%2Fcurrent%2FProjects%2FCLLE-ES%2Fmeth_array.CLLE-ES.tsv', 'download?fn=%2Fcurrent%2FProjects%2FCOAD-US%2Fmeth_array.COAD-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FGBM-US%2Fmeth_array.GBM-US.tsv',
           'download?fn=%2Fcurrent%2FProjects%2FKIRC-US%2Fmeth_array.KIRC-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FKIRP-US%2Fmeth_array.KIRP-US.tsv',
           'download?fn=%2Fcurrent%2FProjects%2FLGG-US%2Fmeth_array.LGG-US.tsv',  'download?fn=%2Fcurrent%2FProjects%2FLIHC-US%2Fmeth_array.LIHC-US.tsv',
          'download?fn=%2Fcurrent%2FProjects%2FOV-US%2Fmeth_array.OV-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FPAEN-AU%2Fmeth_array.PAEN-AU.tsv', 'download?fn=%2Fcurrent%2FProjects%2FPBCA-DE%2Fmeth_array.PBCA-DE.tsv',
          'download?fn=%2Fcurrent%2FProjects%2FPRAD-US%2Fmeth_array.PRAD-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FREAD-US%2Fmeth_array.READ-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FSKCM-US%2Fmeth_array.SKCM-US.tsv',
          'download?fn=%2Fcurrent%2FProjects%2FSTAD-US%2Fmeth_array.STAD-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FTHCA-US%2Fmeth_array.THCA-US.tsv', 'download?fn=%2Fcurrent%2FProjects%2FUCEC-US%2Fmeth_array.UCEC-US.tsv']

files = ['/Users/khandekara2/Documents/methylationProject/01_data/raw_data/ICGC/meth_array.PAEN-AU.tsv']
import csv
for file in files:
    with open(file, 'r') as f, open(file[file.find('.')+1 : file.rfind('.')] + '.tsv', 'w') as csvout:
        reader = csv.reader(f, delimiter='\t')
        csvout = csv.writer(csvout, delimiter='\t')
        rownum = 0
        for row in reader:
            if rownum == 0:
                csvout.writerow(row[4:9])
            elif row[7] in exon_ids:
                csvout.writerow(row[4:9]) 
            rownum += 1
                
               
        


# In[ ]:



