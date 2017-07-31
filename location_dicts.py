import pandas as pd

pbca_dict = {} #key is sample, value is dictionary(location: methylation ratio)
maly_dict = {} #key is sample, value is dictionary(location: methylation ratio)

pbca = pd.read_csv('PBCA_expanded.bed', sep='\t')
maly = pd.read_csv('MALY_expanded.bed', sep='\t')


for sample in pbca[id].unique():
    pbca_dict[sample] = {}
    df = pbca[pbca['id'] == sample]
    for chr, start, stop, ratio in zip(df['chromosome'], df['start'], df['stop'], df['methylation_ratio']):
        pbca[sample][(chr, start, stop)] = ratio

for sample in maly[id].unique():
    pbca_dict[sample] = {}
    df = maly[maly['id'] == sample]
    for chr, start, stop, ratio in zip(df['chromosome'], df['start'], df['stop'], df['methylation_ratio']):
        maly[sample][(chr, start, stop)] = ratio






