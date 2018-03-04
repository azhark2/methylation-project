# create a dictionary that maps the genomic location of a CpG mutation to a list of all the samples that that mutation occurred in

import pickle
import pandas as pd

ssms = ['/data/khandekara2/mutation_data/raw_data/simple_somatic_mutation.open.MALY.tsv', '/data/khandekara2/mutation_data/raw_data/simple_somatic_mutation.open.PBCA.tsv']
common samples = pickle.load(open('common_samples.pickle', 'rb'))

for file in ssms:
	cancer_type = file.split('.')[2]
	mutation_dict = {} #key is chromosomal location, value is list with id's of all samples that the the mutation occurred in 
	df = pd.read_csv(file, sep='\t')
	#process raw mutation file
	df = df[['chromosome', 'chromosome_start', 'chromosome_end', 'submitted_sample_id', 'mutated_from_allele', 'mutated_to_allele', 'total_read_count', 'mutant_allele_read_count',  'consequence_type', 'gene_affected']]
	df['chromosome'] = 'chr' + df['chromosome']
	df.fillna('MISSING', inplace=True)
	df.rename(columns={'submitted_sample_id': 'id'}, inplace=True)
	df1 = df_sample[(df_sample.mutated_from_allele == 'C') & (df_sample.mutated_to_allele == 'T')]
	df2 = df_sample[(df_sample.mutated_from_allele == 'G') & (df_sample.mutated_to_allele == 'A')]
	df3 = pd.concat([df1, df2])
	df3['chromosome_start'].apply(lambda x: x - 1)
	df3 = df3[df3.id in common_samples]
	#add location column
	locations = []
	for chrom, start, stop in zip(list(df3['chromosome']), list(df3['start']), list(df3['stop'])):
    	locations.append((chrom, start, stop))
		df3['location'] = locations
	#populate dictionary
	for loc, sample in zip(df.location, df.id):
		coords = (str(loc[0]), int(loc[1]), int(loc[2]))
		if loc not in mutation_dict:
			mutation_dict[coords] = []
		mutation_dict.append(sample)
			
	#pickle dictionary
	pickle.dump(mutation_dict, open(cancer_type + '_mutation_dict.pickle', 'wb'))
		
		

