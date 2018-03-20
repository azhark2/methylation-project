#filter original colorectal cancer mutation calls for only samples that were classified as POLE-mutants

import pandas as pd

df = pd.read_csv('Poulos_NAR_2017_colorectal_mutations.tsv', sep='\t')
df['reference'] = df['reference'].apply(lambda x: x.upper())
df.id.apply(str)
POLE_samples = set(['TCGA-AZ-4315', 'TCGA-F5-6814', 'TCGA-CA-6717', 'TCGA-AA-3977', 'TCGA-CA-6718', 'TCGA-AA-3555', 'TCGA-A6-6141'])


if df['mutated_from'].equals(df['reference']):
	# pole = df[df['id'] in POLE_samples]
	pole = df[df['id'].str.contains('TCGA-AZ-4315') | df['id'].str.contains('TCGA-F5-6814') | df['id'].str.contains('TCGA-CA-6717') | df['id'].str.contains('TCGA-AA-3977') | df['id'].str.contains('TCGA-CA-6718') | df['id'].str.contains('TCGA-AA-3555') | df['id'].str.contains('TCGA-A6-6141')]
	pole.to_csv('POLE_mutants_only.tsv', sep='\t', index=False, header=False)
