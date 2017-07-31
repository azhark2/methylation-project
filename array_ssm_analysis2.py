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
    mutations = []
    for chromosome, position in zip(chromosomes, positions):
        mutations.append((chromosome, position))
    #now look for overlaps
    overlaps = list(set(mutations) & set(cpg_sites))
    #store overlaps
    cancertypeToOverlaps[cancer_type] = overlaps




ssms = ['download?fn=%2Frelease_23%2FProjects%2FPAEN-AU%2Fsimple_somatic_mutation.open.PAEN-AU.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBLCA-US%2Fsimple_somatic_mutation.open.BLCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBRCA-US%2Fsimple_somatic_mutation.open.BRCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCESC-US%2Fsimple_somatic_mutation.open.CESC-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCOAD-US%2Fsimple_somatic_mutation.open.COAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FGBM-US%2Fsimple_somatic_mutation.open.GBM-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FKIRC-US%2Fsimple_somatic_mutation.open.KIRC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FKIRP-US%2Fsimple_somatic_mutation.open.KIRP-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FLGG-US%2Fsimple_somatic_mutation.open.LGG-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FLIHC-US%2Fsimple_somatic_mutation.open.LIHC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FOV-US%2Fsimple_somatic_mutation.open.OV-US.tsv',
         'download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv', 'download?fn=%2Frelease_23%2FProjects%2FREAD-US%2Fsimple_somatic_mutation.open.READ-US.tsv'
        'download?fn=%2Frelease_23%2FProjects%2FSKCM-US%2Fsimple_somatic_mutation.open.SKCM-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FSTAD-US%2Fsimple_somatic_mutation.open.STAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FTHCA-US%2Fsimple_somatic_mutation.open.THCA-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FUCEC-US%2Fsimple_somatic_mutation.open.UCEC-US.tsv']




cancertypeToOverlaps = {} #dictionary that maps cancer type to list of overlaps
sample_ToOverlaps = {}
ssm_data = {}

manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
df2 = pd.read_csv(manifest_location, header=7)

#create list of tuples of all assayed cpg site locations
cpg_sites = []
for c, p in zip(list(df2['CHR']), list(df2['MAPINFO'])):
    cpg_sites.append((c,p))

for s in ssms:
    cancer_type = s.split('.')[2]
    ssm = pd.read_csv(s, sep='\t')
    ssm = ssm[
        (ssm.reference_genome_allele == 'C') & (ssm.mutated_to_allele == 'T') | (ssm.reference_genome_allele == 'G') & (
        ssm.mutated_to_allele == 'A')]
    totalCtoT = ssm.shape[0]
    ssm = ssm[ssm.consequence_type == 'missense_variant']
    cdsCtoT = ssm.shape[0]
    ssm_data[cancer_type] = (totalCtoT, cdsCtoT)
    num_sample_mut = ssm.submitted_sample_id.nunique()
    calculate_overlaps()




pickle.dump(cancertypeToOverlaps, open("cancertypeToOverlaps", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
# pickle.dump(sample_ToOverlaps, open("sample_ToOverlaps.pickle", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
print (cancertypeToOverlaps.items())