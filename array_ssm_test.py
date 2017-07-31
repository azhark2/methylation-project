
# coding: utf-8

# In[13]:

from __future__ import division
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
import pickle

def get_avgBetaValue():
    groupby = array['methylation_value'].groupby(array['probe_id'])
    df = groupby.agg([np.mean, np.std])
    df.reset_index(inplace=True) #creates dataframe with probe_id, mean beta value, and std deviation as columns
    cancertypeToAvgBeta[cancer_type] = {} #dictionary with cancer type as key and dictionary that maps cg to average beta as value
    for cg, avg, std in zip(list(df['probe_id']), list(df['mean']), list(df['std'])):
        cancertypeToAvgBeta[cancer_type][cg] = (avg, std)


#calculate # of CpG sites that are both mutated and methylated across samples and store information in dictionary
def calculate_overlaps():
    chromosomes = ssm['chromosome']
    positions = ssm['chromosome_start']
    types = ssm['consequence_type']
    mutations = []
    cpg_sites = []
    for chromosome, position in zip(chromosomes, positions):
        mutations.append((chromosome, position))
    for cg_id in array['probe_id']:
        cpg_sites.append(probeToLocation[cg_id])
    #now look for overlaps
    overlaps = list(set(mutations) & set(cpg_sites))
    #store overlaps
    cancertypeToOverlaps[cancer_type] = overlaps

def calculate_sampleOverlaps():
    array_samples = list(array.submitted_sample_id.unique())
    ssm_samples = list(ssm.submitted_sample_id.unique())
    sample_overlaps = list(set(array_samples) & set(ssm_samples))
    for sample in sample_overlaps:
        sample_array = array[array.submitted_sample_id == sample]
        sample_ssms = ssm[ssm.submitted_sample_id == sample]
        chromosomes = sample_ssms['chromosome']
        positions = sample_ssms['chromosome_start']
        mutations = []
        cpg_sites = []
        for chromosome, position in zip(chromosomes, positions):
            mutations.append((chromosome, position))
        for cg_id in (sample_array.probe_id):
            cpg_sites.append(probeToLocation[cg_id])
        overlaps = set(set(mutations) & set(cpg_sites))
        sample_ToOverlaps[sample] = overlaps
        methylated_overlaps = []
        #filter array dataframe for only methylated sites
        df = sample_array[sample_array['methylation_value'] > 0.7]
        for cg_id, beta_value in zip(list(df['probe_id']), list(df['methylation_value'])):
            if probeToLocation[cg_id] in overlaps:
                methylated_overlaps.append(probeToLocation[cg_id])
        sample_ToMethylatedOverlaps[sample] = methylated_overlaps
#main  
manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
manifest = pd.read_csv(manifest_location, header=7)
probeToLocation = {}
for id, chr, coord in zip(list(manifest.IlmnID.dropna()), list(manifest['CHR'].dropna()), list(manifest['MAPINFO'].dropna())):
    probeToLocation[id] = (chr, int(coord))
arrays = ['PAEN-AU.tsv']
ssms = ['/Users/khandekara2/Documents/methylationProject/01_data/raw_data/ICGC/simple_somatic_mutation.open.PAEN-AU.tsv']
cancertypeToOverlaps = {} #dictionary that maps cancer type to list of overlaps
sample_ToOverlaps = {} #dictionary that maps sample to list of overlaps
sample_ToMethylatedOverlaps = {}
cancertypeToAvgBeta = {}
cancertypeToMethylatedOverlaps = {}
ssm_data = {} #dictionary that maps cancer type to #of total CtoT/GtoA mutations,
for a, s in zip(arrays, ssms):
    cancer_type = a.split('.')[0]
    array = pd.read_csv(a, sep='\t')
    ssm = pd.read_csv(s, sep='\t')
    ssm = ssm[(ssm.reference_genome_allele == 'C') & (ssm.mutated_to_allele == 'T') | (ssm.reference_genome_allele == 'G') & (
        ssm.mutated_to_allele == 'A')]
    totalCtoT = ssm.shape[0]
    ssm = ssm[ssm.consequence_type == 'missense_variant']
    cdsCtoT = ssm.shape[0]
    ssm_data[cancer_type] = (totalCtoT, cdsCtoT)
    num_sample_mut = ssm.submitted_sample_id.nunique()
    num_sample_arr = array.submitted_sample_id.nunique()
    get_avgBetaValue()
    calculate_overlaps()
    calculate_sampleOverlaps()


# In[14]:

print (sample_ToOverlaps.items())


# In[15]:

print (cancertypeToOverlaps.items())


# In[10]:

print (cancertypeToAvgBeta.items())


# In[19]:

get_ipython().system('pip install cPickle')


# In[4]:

import pickle as cPickle
import numpy

cancertypeToAvgBeta = cPickle.load(open("cancertypeToAvgBeta.pickle", "rb"))
sample_ToOverlaps = cPickle.load(open("sample_ToOverlaps.pickle", "rb"))
cancertypeToOverlaps = cPickle.load(open("cancertypeToOverlaps.pickle", "rb"))



# In[5]:

print (sample_ToOverlaps.items())


# In[6]:

print (cancertypeToOverlaps)


# In[9]:

for cancer, overlaps in dict.items():
    if len(dict[cancer]) != 0:
        print (cancer, overlaps)


# In[48]:

for cancer, num in cancertypeToOverlaps.items():
    print (cancer + len(num))
    


# In[ ]:

from __future__ import division
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
import pickle

sys.argsv = []

# def get_cancerType(a):
#     mystr = a
#     char1 = '.'
#     char2 = '.'
#     return mystr[mystr.find(char1)+1 : mystr.rfind(char2)]


def get_avgBetaValue():
    groupby = array['methylation_value'].groupby(array['probe_id'])
    df = groupby.agg([np.mean, np.std])
    df.reset_index(inplace=True) #creates dataframe with probe_id, mean beta value, and std deviation as columns
    cancertypeToAvgBeta[cancer_type] = {} #dictionary with cancer type as key and dictionary that maps cg to average beta as value
    for cg, avg, std in zip(list(df['probe_id']), list(df['mean']), list(df['std'])):
        cancertypeToAvgBeta[cancer_type][cg] = (avg, std)


#calculate # of CpG sites that are both mutated and methylated across samples and store information in dictionary
def calculate_overlaps():
    chromosomes = ssm['chromosome']
    positions = ssm['chromosome_start']
    types = ssm['consequence_type']
    mutations = []
    cpg_sites = []
    for chromosome, position in zip(chromosomes, positions):
        mutations.append((chromosome, position))
    for cg_id in array['probe_id']:
        cpg_sites.append(probeToLocation[cg_id])
    #now look for overlaps
    overlaps = list(set(mutations) & set(cpg_sites))
    #store overlaps
    cancertypeToOverlaps[cancer_type] = overlaps

def calculate_sampleOverlaps():
    array_samples = list(array.submitted_sample_id.unique())
    ssm_samples = list(ssm.submitted_sample_id.unique())
    sample_overlaps = list(set(array_samples) & set(ssm_samples))
    for sample in sample_overlaps:
        sample_array = array[array.submitted_sample_id == sample]
        sample_ssms = ssm[ssm.submitted_sample_id == sample]
        chromosomes = sample_ssms['chromosome']
        positions = sample_ssms['chromosome_start']
        mutations = []
        cpg_sites = []
        for chromosome, position in zip(chromosomes, positions):
            mutations.append((chromosome, position))
        for cg_id in (sample_array.probe_id):
            cpg_sites.append(probeToLocation[cg_id])
        overlaps = set(set(mutations) & set(cpg_sites))
        sample_ToOverlaps[sample] = overlaps
        methylated_overlaps = []
        #filter array dataframe for only methylated sites
        df = sample_array[sample_array['methylation_value'] > 0.7]
        for cg_id, beta_value in zip(list(df['probe_id']), list(df['methylation_value'])):
            if probeToLocation[cg_id] in overlaps:
                methylated_overlaps.append(probeToLocation[cg_id])
        sample_ToMethylatedOverlaps[sample] = methylated_overlaps



ssms = ['download?fn=%2Frelease_23%2FProjects%2FPAEN-AU%2Fsimple_somatic_mutation.open.PAEN-AU.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBLCA-US%2Fsimple_somatic_mutation.open.BLCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBRCA-US%2Fsimple_somatic_mutation.open.BRCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCESC-US%2Fsimple_somatic_mutation.open.CESC-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCOAD-US%2Fsimple_somatic_mutation.open.COAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FGBM-US%2Fsimple_somatic_mutation.open.GBM-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FKIRC-US%2Fsimple_somatic_mutation.open.KIRC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FKIRP-US%2Fsimple_somatic_mutation.open.KIRP-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FLGG-US%2Fsimple_somatic_mutation.open.LGG-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FLIHC-US%2Fsimple_somatic_mutation.open.LIHC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FOV-US%2Fsimple_somatic_mutation.open.OV-US.tsv',
         'download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv', 'download?fn=%2Frelease_23%2FProjects%2FREAD-US%2Fsimple_somatic_mutation.open.READ-US.tsv'
        'download?fn=%2Frelease_23%2FProjects%2FSKCM-US%2Fsimple_somatic_mutation.open.SKCM-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FSTAD-US%2Fsimple_somatic_mutation.open.STAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FTHCA-US%2Fsimple_somatic_mutation.open.THCA-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FUCEC-US%2Fsimple_somatic_mutation.open.UCEC-US.tsv']




cancertypeToOverlaps = {} #dictionary that maps cancer type to list of overlaps
sample_ToOverlaps = {}


manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
df2 = pd.read_csv(manifest_location, header=7)

#create list of tuples of all assayed cpg site locations
cpg_sites = []
for c, p in zip(list(df2['CHR']), list(df2['MAPINFO'])):
    cpg_sites.append((c,p))

for s in ssms:
    cancer_type = a.split('.')[2]
    ssm = pd.read_csv(s, sep='\t')
    ssm = ssm[
        (ssm.reference_genome_allele == 'C') & (ssm.mutated_to_allele == 'T') | (ssm.reference_genome_allele == 'G') & (
        ssm.mutated_to_allele == 'A')]
    totalCtoT = ssm.shape[0]
    ssm = ssm[ssm.consequence_type == 'missense_variant']
    cdsCtoT = ssm.shape[0]
    ssm_data[cancer_type] = (totalCtoT, cdsCtoT)
    num_sample_mut = ssm.submitted_sample_id.nunique()
    num_sample_arr = array.submitted_sample_id.nunique()
    calculate_overlaps()
    calculate_sampleOverlaps()



pickle.dump(cancertypeToOverlaps, open("cancertypeToOverlaps", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(sample_ToOverlaps, open("sample_ToOverlaps.pickle", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

