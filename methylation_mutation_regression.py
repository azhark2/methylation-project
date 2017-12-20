#simple linear regression to model the relationship between methylation ratio and mutation frequency in various cancers
#goal is to compare the models derived from coding regions vs. all other regions
#normal WGBS data from various tissues is used in conjunction with mutation data of cancers originating in same tissue 

import pybedtools

methylation_files = ['liver', 'esophagus', 'breast', 'lung', 'ovary', 'pancreas', 'thymus']
mutation_files = ['LIHC-US', 'ESAD-UK', 'LUAD-US', 'OV', 'PAEN-AU', ]
cds_files = [] #normal WGBS intersected with cds

# cds = pybedtools.BedTool(cds.bed)
for meth, mut, cds_meth in zip(methylation_files, mutation_files, cds_files):
	a = pybedtools.BedTool(meth)
	b = pybedtools.BedTool(cds_meth)
	a.intersect(mut).saveas('' + 'methylation_mutation.bed')
	b.intersect(mut).saveas('' + 'methylation_mutation_cds.bed')
		
	
	