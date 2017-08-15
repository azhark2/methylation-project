import csv
#concatenates individual bed files into one large bed file
# bed_files = ['chr1_cds_cpg.bed', 'chr2_cds_cpg.bed', 'chr3_cds_cpg.bed', 'chr4_cds_cpg.bed', 'chr5_cds_cpg.bed', 'chr6_cds_cpg.bed', 'chr7_cds_cpg.bed'
#             'chr8_cds_cpg.bed', 'chr9_cds_cpg.bed', 'chr10_cds_cpg.bed', 'chr11_cds_cpg.bed', 'chr12_cds_cpg.bed', 'chr13_cds_cpg.bed', 'chr14_cds_cpg.bed'
#             'chr15_cds_cpg.bed', 'chr16_cds_cpg.bed', 'chr17_cds_cpg.bed', 'chr18_cds_cpg.bed', 'chr19_cds_cpg.bed', 'chr20_cds_cpg.bed',
#             'chr21_cds_cpg.bed', 'chr22_cds_cpg.bed']


# prefix = 'CpG_hg19_chr'
# suffix = '.bed'
# for i in range(1, 23):
#     bed_files.append(prefix + str(i) + suffix)

bed_files = ['meth_seq_CLLE_cds.bed', 'meth_seq_MALY_cds.bed', 'meth_seq_PBCA_cds.bed']
with open('meth_seq_cds', 'w') as csvout:
    writer = csv.writer(csvout, delimiter='\t')
    for file in bed_files:
        with open(file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                writer.writerow(row)


