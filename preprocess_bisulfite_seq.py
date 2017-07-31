
#preprocesses whole genome bisulfite sequencing files from ICGC by filtering for only CpG's that are contained in a cds(protein coding) region
import csv
import pickle
# list_of_coding_cpgs = set([]) #coordinates of coding cpg's
# bed_files = ['chr1_cds_cpg.bed', 'chr2_cds_cpg.bed', 'chr3_cds_cpg.bed', 'chr4_cds_cpg.bed', 'chr5_cds_cpg.bed', 'chr6_cds_cpg.bed', 'chr7_cds_cpg.bed',
#             'chr8_cds_cpg.bed', 'chr9_cds_cpg.bed', 'chr10_cds_cpg.bed', 'chr11_cds_cpg.bed', 'chr12_cds_cpg.bed', 'chr13_cds_cpg.bed', 'chr14_cds_cpg.bed',
#             'chr15_cds_cpg.bed', 'chr16_cds_cpg.bed', 'chr17_cds_cpg.bed', 'chr18_cds_cpg.bed', 'chr19_cds_cpg.bed', 'chr20_cds_cpg.bed',
#             'chr21_cds_cpg.bed', 'chr22_cds_cpg.bed']


# with open('all_cds_cpgs', 'r') as f:
#     reader = csv.reader(f, delimiter='\t')
#     for row in reader:
#         list_of_coding_cpgs.add((row[0], row[1]))
#         list_of_coding_cpgs.add((row[0], row[2]))
# pickle.dump(list_of_coding_cpgs, open('list_of_coding_cpgs.pickle', 'wb'))

#transform bisulfite sequencing files into bed format
biseq_files = ['meth_seq_MALY.tsv', 'meth_seq_CLLE.tsv', 'meth_seq_PBCA.tsv']
for file in biseq_files:
    with open(file, 'r') as f, open(file.split('.')[0] + '.bed', 'w') as csvout:
        reader = csv.reader(f, delimiter='\t')
        csvout = csv.writer(csvout, delimiter='\t')
        rownum = 0
        for row in reader:
            if rownum > 0:
                csvout.writerow(['chr' + str(row[6]), int(row[7]), int(row[8]), row[4], row[11], row[12]]) #chrom, start, stop, submitted_sample_id, methylation_ratio
            rownum += 1

