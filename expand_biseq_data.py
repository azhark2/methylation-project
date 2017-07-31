#preprocess bisulfite sequencing data so that methylation information from one cytosine is "expanded" to neighboring cytosine
import csv

biseq_files = ['meth_seq_CLLE_cds.bed', 'meth_seq_MALY_cds.bed', 'meth_seq_PBCA_cds.bed']
for file in biseq_files:
    with open(file, 'r') as f, open(file.split('_')[2] + 'expanded.bed', 'w') as csvout:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(csvout, delimiter='\t')
        rownum = 0
        for row in reader:
            writer.writerow([row[0], row[1], row[2], row[3], row[4]]) #chrom, start, stop, id, methylation value
            writer.writerow([row[0], int(row[1]) + 1, int(row[2]) + 1, row[3], row[4]]) #chrom, start, stop, id, methylation value


