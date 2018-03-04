#create dictionary that stores methylation ratio for every site in every sample in PBCA and MALY
import pickle
import pandas as pd


pbca_dict = {} #key is sample, value is dictionary(location: methylation ratio)
maly_dict = {} #key is sample, value is dictionary(location: methylation ratio)

df = pd.read_csv('PBCA_cds_expanded.tsv', sep='\t')
for sample in df['id'].unique():
    pbca_dict[sample] = {}
    df1 = df[df['id'] == sample]
    for chrom, start, stop, ratio in zip(df1['chromosome'], df1['start'], df1['stop'], df1['methylation_ratio']):
        pbca_dict[sample][(str(chrom), int(start), int(stop))] = ratio

df = pd.read_csv('MALY_cds_expanded.tsv', sep='\t')
for sample in df['id'].unique():
    maly_dict[sample] = {}
    df1 = df[df['id'] == sample]
    for chrom, start, stop, ratio in zip(df1['chromosome'], df1['start'], df1['stop'], df1['methylation_ratio']):
        maly_dict[sample][(str(chrom), int(start), int(stop))] = ratio


pickle.dump(pbca_dict, open('pbca_dict.pickle,', 'wb'))
pickle.dump(maly_dict, open('pbca_dict.pickle,', 'wb'))

