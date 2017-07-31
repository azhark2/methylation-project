
#create bed files of all somatic point mutations in cds region
import pandas as pd
import pybedtools

ssms = ['simple_somatic_mutation.open.CLLE-ES.tsv', 'simple_somatic_mutation.open.MALY-DE.tsv', 'simple_somatic_mutation.open.PBCA-DE.tsv']

for file in ssms:
    df = pd.read_csv(file, sep='\t')
    df = df[['chromosome', 'chromosome_start', 'chromosome_end', 'submitted_sample_id', 'mutated_from_allele',
             'mutated_to_allele', 'consequence_type', 'gene_affected']]
    df.fillna('MISSING', inplace=True)
    df.to_csv(file[0:-3] + '.raw.bed', sep='\t', index=False, header=False)
    df.to_csv(file[0:-3] + '.raw.bed', sep='\t', index=False, header=False)

    a = pybedtools.BedTool(file)
    a.intersect('cds.bed').saveas(file[0:-3] + 'bed')



