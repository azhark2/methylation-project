#write out ICGC somatic mutation files as bed files
import pandas as pd

all_ssms = ['download?fn=%2Frelease_23%2FProjects%2FPAEN-AU%2Fsimple_somatic_mutation.open.PAEN-AU.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBLCA-US%2Fsimple_somatic_mutation.open.BLCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FBRCA-US%2Fsimple_somatic_mutation.open.BRCA-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCESC-US%2Fsimple_somatic_mutation.open.CESC-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCOAD-US%2Fsimple_somatic_mutation.open.COAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FGBM-US%2Fsimple_somatic_mutation.open.GBM-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FKIRC-US%2Fsimple_somatic_mutation.open.KIRC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FKIRP-US%2Fsimple_somatic_mutation.open.KIRP-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FLGG-US%2Fsimple_somatic_mutation.open.LGG-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FLIHC-US%2Fsimple_somatic_mutation.open.LIHC-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FOV-US%2Fsimple_somatic_mutation.open.OV-US.tsv',
         'download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv', 'download?fn=%2Frelease_23%2FProjects%2FREAD-US%2Fsimple_somatic_mutation.open.READ-US.tsv'
        'download?fn=%2Frelease_23%2FProjects%2FSKCM-US%2Fsimple_somatic_mutation.open.SKCM-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FSTAD-US%2Fsimple_somatic_mutation.open.STAD-US.tsv', 'download?fn=%2Frelease_23%2FProjects%2FTHCA-US%2Fsimple_somatic_mutation.open.THCA-US.tsv',
        'download?fn=%2Frelease_23%2FProjects%2FUCEC-US%2Fsimple_somatic_mutation.open.UCEC-US.tsv']

# ssms = ['simple_somatic_mutation.open.CLLE-ES.tsv', 'simple_somatic_mutation.open.MALY-DE.tsv', 'simple_somatic_mutation.open.PBCA-DE.tsv']

for file in all_ssms:
    df = pd.read_csv(file, sep='\t')
    df1 = df[(df.mutated_from_allele == 'C') & (df.mutated_to_allele == 'T')]
    df2 = df[(df.mutated_from_allele == 'G') & (df.mutated_to_allele == 'A')]
    df = pd.concat([df1, df2])
    cols = list(df.columns)
    df = df[['chromosome', 'chromosome_start', 'chromosome_end', 'submitted_sample_id', 'mutated_from_allele', 'mutated_to_allele',  'consequence_type', 'gene_affected']]
    # df['chromosome'] = 'chr' + df['chromosome']
    df.fillna('MISSING', inplace=True)
    df.rename(columns={'submitted_sample_id': 'id'}, inplace=True)
    df.to_csv(file.split('.')[1] + '_mutation.bed.all', sep='\t', index=False)


