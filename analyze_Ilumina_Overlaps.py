
# coding: utf-8

# In[1]:

import pickle as cPickle
import numpy

cancertypeToAvgBeta = cPickle.load(open("cancertypeToAvgBeta.pickle", "rb"))
cancertypeToOverlaps = cPickle.load(open("overlaps.pickle", "rb"))


# In[13]:

print (cancertypeToAvgBeta.items()[0])


# In[ ]:




# In[12]:

import pandas as pd
manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
manifest = pd.read_csv(manifest_location, header=7)
probeToLocation = {}
for id, chr, coord in zip(list(manifest.IlmnID.dropna()), list(manifest['CHR'].dropna()), list(manifest['MAPINFO'].dropna())):
    probeToLocation[id] = (chr, int(coord))


# In[14]:

import csv
with open('overlaps.csv', 'w') as csvout:
    csvout = csv.writer(csvout, delimiter='\t')
    for cancer, overlap in cancertypeToOverlaps.items():
        if len(cancertypeToOverlaps[cancer]) != 0:
            print (cancertypeToAvgBeta[cancer][overlap][0]])
#             csvout.writerow([cancer, overlap, cancertypeToAvgBeta[cancer][overlap][0]])


# In[ ]:



