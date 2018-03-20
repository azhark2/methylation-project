#filter all files for only C>T or G>A mutations(anywhere in genome)
for file in os.listdir("/data/khandekara2/mutation_data/raw_data"):
    if file.endswith(".tsv"):
    	df = pd.read_csv(file, sep='\t')
    	df1 = df[(df.mutated_from_allele == 'C') & (df.mutated_to_allele == 'T')]
    	df2 = df[(df.mutated_from_allele == 'G') & (df.mutated_to_allele == 'A')]
    	df = pd.concat([df1, df2])
    	cols = list(df.columns)
    	df = df[['chromosome', 'chromosome_start', 'chromosome_end', 'submitted_sample_id', 'mutated_from_allele', 'mutated_to_allele', 'total_read_count', 'mutant_allele_read_count',  'consequence_type', 'gene_affected']]
    	df.chromosome = df.chromosome.astype(str)
        df['chromosome'] = 'chr' + df['chromosome']
    	df.fillna('MISSING', inplace=True)
    	df.rename(columns={'submitted_sample_id': 'id'}, inplace=True)
    	df.to_csv(file.split('.')[2] + '_mutation.bed.all', sep='\t', index=False, header=False)
        remove_duplicates.remove(file.split('.')[2] + '_mutation.bed.all')
