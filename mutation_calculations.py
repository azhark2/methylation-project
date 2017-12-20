
#script to calculate total # of mutations, # of mutations in coding regions, and # of mutations at CpG sites in coding regions
import pandas as pd
import csv
import pybedtools
import os
import remove_duplicates

# all_ssms = ['download?fn=%2Frelease_23%2FProjects%2FPAEN-AU%2Fsimple_somatic_mutation.open.PAEN-AU.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBLCA-US%2Fsimple_somatic_mutation.open.BLCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBRCA$
#         'download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCOAD-US%2Fsimple_somatic_mutation.open.COAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FGBM-US%2$
#         'download?fn=%2Frelease_23%2FProjects%2FKIRC-US%2Fsimple_somatic_mutation.open.KIRC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FKIRP-US%2Fsimple_somatic_mutation.open.KIRP-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FLGG-US%2$
#         'download?fn=%2Frelease_23%2FProjects%2FLIHC-US%2Fsimple_somatic_mutation.open.LIHC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FOV-US%2Fsimple_somatic_mutation.open.OV-US.tsv',
#          'download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv', 'download?fn=%2Frelease_23%2FProjects%2FREAD-US%2Fsimple_somatic_mutation.open.READ-US.tsv', 'simple_somatic_mutation.open.CLLE-ES.tsv', 'si$
#         'download?fn=%2Frelease_23%2FProjects%2FSKCM-US%2Fsimple_somatic_mutation.open.SKCM-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FSTAD-US%2Fsimple_somatic_mutation.open.STAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FTHCA-US%$
#         'download?fn=%2Frelease_23%2FProjects%2FUCEC-US%2Fsimple_somatic_mutation.open.UCEC-US.tsv']

ssms = ['/data/khandekara2/mutation_data/raw_data/simple_somatic_mutation.open.CLLE.tsv', '/data/khandekara2/mutation_data/raw_data/simple_somatic_mutation.open.MALY.tsv', '/data/khandekara2/mutation_data/raw_data/simple_somatic_mutation.open.PBCA.tsv']

#filter all files for only C>T or G>A mutations(anywhere in genome)
for file in ssms:
	with open(file.split('.')[2] + '_mutation_calculations.tsv', 'w') as csvout:
		writer = csv.writer(csvout, delimiter='\t')
		cancer_type = file.split('.')[2]
		#process raw mutation file
    	df = pd.read_csv(file, sep='\t')
    	df = df[['chromosome', 'chromosome_start', 'chromosome_end', 'submitted_sample_id', 'mutated_from_allele', 'mutated_to_allele', 'total_read_count', 'mutant_allele_read_count',  'consequence_type', 'gene_affected']]
    	df['chromosome'] = 'chr' + df['chromosome']
    	df.fillna('MISSING', inplace=True)
    	df.rename(columns={'submitted_sample_id': 'id'}, inplace=True)
    	#calculate total number of mutations for each sample
    	count = 0
    	for sample in df['id'].unique():
    		count += 1
    		num_total = 0
    		num_CpG = 0
    		num_CpG_coding = 0	
    		df_sample = df[df['id'] == sample]
    	    #calculate total number of mutations
    	    df_sample.to_csv(sample + '_all_mutations.bed', sep='\t', index=False, header=False)
    	    #remove duplicates and rename file with no duplicates to original file
    	    remove_duplicates.remove(sample + '_all_mutations.bed')
    	    os.rename(sample + '_all_mutations.bed' + '.noDuplicates', sample + '_all_mutations.bed')
    	    num_total = sum(1 for line in open(sample + '_all_mutations.bed'))
    	    #calculate total number of mutations in coding regions
    	    a = pybedtools.BedTool(sample + '_all_mutations.bed')
    	    a.intersect('cds.bed').saveas(sample + '_coding_mutations.bed')
    	    remove_duplicates.remove(sample + '_coding_mutations.bed')
    	    os.rename(sample + '_coding_mutations.bed' + '.noDuplicates', sample + '_coding_mutations.bed')
    	    num_coding_mutations = sum(1 for line in open(sample + '_coding_mutations.bed'))
    		#calculate total number of CpG mutations per sample
    		df1 = df_sample[(df_sample.mutated_from_allele == 'C') & (df_sample.mutated_to_allele == 'T')]
    		df2 = df_sample[(df_sample.mutated_from_allele == 'G') & (df_sample.mutated_to_allele == 'A')]
    		df3 = pd.concat([df1, df2])    		 
    		df3.to_csv(sample + '_all_CpG_mutations.bed', sep='\t', index=False, header=False)
    		a = pybedtools.BedTool(sample + '_all_CpG_mutations.bed')
    		a.intersect('all_cpgs_shifted.bed').saveas(sample + '_all_CpG_mutations.bed'')
    		remove_duplicates.remove(sample + '_all_CpG_mutations.bed')
    	    os.rename(sample + '_all_CpG_mutations.bed' + '.noDuplicates', sample + '_all_CpG_mutations.bed')
    		num_CpG = sum(1 for line in open(sample + '_all_CpG_mutations.bed')) #total number of CpG mutations per sample
    		#calculate # of CpG mutations in coding regions only
    		a = pybedtools.BedTool(sample + '_all_CpG_mutations.bed')
    		a.intersect('cds.bed').saveas(sample + '_coding_CpG_mutations') #file of all CpG mutations in coding region for that sample
    		remove_duplicates.remove(sample + '_coding_CpG_mutations.bed')
    	    os.rename(sample + '_mutation.bed.cds' + '.noDuplicates', sample + '_coding_CpG_mutations.bed')
    		num_CpG_coding = sum(1 for line in open(sample + '_coding_CpG_mutations.bed')) #number of CpG mutations in coding region
    		writer.writerow([cancer_type, sample, num_total, num_CpG, num_coding_mutations, num_CpG_coding]) #submitted_sample_id, total # of mutations, #of coding mutations, # of CpG mutations, # of CpG mutations in coding regions
    		#for testing purposes
    		if count > 3:
    			break
    		
    