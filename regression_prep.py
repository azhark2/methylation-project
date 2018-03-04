#simple linear regression to model the relationship between methylation ratio and mutation frequency in various cancers
#goal is to determine how much methylation expalins the variation in coding and non-coding regions compare the models derived from coding regions vs. all other regions
#normal WGBS data from various tissues is used in conjunction with mutation data of cancers originating in same tissue

import pybedtools
# matching of WGBS data from normal tissue to pooled mutation in the same tissue
methylation_files = ['GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig.bed.non_cds',  'GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.wig.bed.non_cds', 'GSM983647_UCSD.Lung.Bisulfite-Seq.STL002.wig.bed.non_cds', 'GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.bed.non_cds',  'GSM983651_UCSD.Pancreas.Bisulfite-Seq.STL003.wig.bed.non_cds', 'GSM1010984_UCSD.Gastric.Bisulfite-Seq.STL003.wig.bed.non_cds', 'GSM1010979_UCSD.Thymus.Bisulfite-Seq.STL001.wig.bed.non_cds', 'normal_colon_WGBS_processed_expanded.bed.non_cds', 'GSM916052_BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig.bed.non_cds']
mutation_files = ['LIVER-COMBINED_mutation.bed.all.noDuplicates', 'BRCA-US_mutation.bed.all.noDuplicates', 'LUNG-COMBINED_mutation.bed.all.noDuplicates', 'OV-US_mutation.bed.all.noDuplicates', 'PANCREAS-COMBINED_mutation.bed.all.noDuplicates', 'GASTRIC-COMBINED_mutation.bed.all.noDuplicates', 'MALY-DE_mutation.bed.all.noDuplicates', 'COAD-US_mutation.bed.all.noDuplicates', 'CLLE-ES_mutation.bed.all.noDuplicates', ]
cds_files = ['GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig.bed.cds', 'GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.wig.bed.cds', 'GSM983647_UCSD.Lung.Bisulfite-Seq.STL002.wig.bed.cds', 'GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.bed.cds', 'GSM983651_UCSD.Pancreas.Bisulfite-Seq.STL003.wig.bed.cds', 'GSM1010984_UCSD.Gastric.Bisulfite-Seq.STL003.wig.bed.cds', 'GSM1010979_UCSD.Thymus.Bisulfite-Seq.STL001.wig.bed.cds', 'normal_colon_WGBS_processed_expanded.bed.cds', 'GSM916052_BI.Mobilized_CD34_Primary_Cells.Bisulfite-Seq.RO_01549.wig.bed.cds'] #normal WGBS intersected with cds

cds = pybedtools.BedTool('/data/khandekara2/bed_CpGs/cds.bed')
for meth, mut, cds_meth in zip(methylation_files, mutation_files, cds_files):
	a = pybedtools.BedTool(meth)
	b = pybedtools.BedTool(cds_meth)
	a.intersect(mut).saveas(mut.split('_')[0] + 'methylation_mutation.bed') #non-cds
	b.intersect(mut).saveas(mut.split('_')[0]  + 'methylation_mutation_cds.bed')

# print(len(methylation_files) == len(mutation_files) == len(cds_files))