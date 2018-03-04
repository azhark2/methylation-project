import pandas as pd

file = 'copy_number_somatic_mutation.BOCA-FR.tsv'
df = pd.read_csv(file, sep='\t')
df = df[['chromosome', 'chromosome_start', 'chromosome_end', 'submitted_sample_id', 'mutation_type', 'copy_number', 'segment_mean',  'segment_median', 'gene_affected']]
df.fillna('MISSING', inplace=True)
df.rename(columns={'submitted_sample_id': 'id'}, inplace=True)
df.to_csv('copy_number_BOCA-FR.tsv', sep='\t', index=False)