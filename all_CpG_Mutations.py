#create a bed file of all CpG mutations
import pandas as pd
import pybedtools
import csv

ssm_files = ['CLLE-ES_mutation.tsv', 'PBCA-DE_mutation.tsv', 'MALY-DE_mutation.tsv']
wgbs_files = ['CLLE_expanded.tsv', 'PBCA_expanded.tsv', 'MALY_expanded.tsv']
frames = []
for wgbs_file, ssm_file in zip(wgbs_files, ssm_files):
    biseq = pd.read_csv(wgbs_file, sep='\t')
    mutation = pd.read_csv(ssm_file, sep='\t')

    biseq_samples = set(biseq['id'].unique())
    mutation_samples = set(mutation['id'].unique())

    #generate list of matching samples in the two types of data
    common_samples = biseq_samples & mutation_samples

    for sample in common_samples:
        sub_mutation = mutation[mutation['id'] == sample]
        frames.append(sub_mutation)

    df = pd.concat(frames)
    df.to_csv('all_CpG_mutations.bed', sep='\t', index=False, header=False)

a = pybedtools.BedTool('all_CpG_mutations.bed')
a.intersect('all_cds_cpgs_shifted_2.bed').saveas('all_cds_CpG_mutations.bed')

#remove duplicates
duplicates = set([])
count = 0

with open('all_cds_CpG_mutations.bed', 'r') as f, open('all_cds_CpG_mutations.noDuplicates', 'w') as csvout:
    reader = csv.reader(f, delimiter='\t')
    writer = csv.writer(csvout, delimiter='\t')
    for row in reader:
        info = (row[0], row[1], row[2], row[3])
        if info not in duplicates:
            duplicates.add(info)
            writer.writerow(row)
        else:
            count += 1




