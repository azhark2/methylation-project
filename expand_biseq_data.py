#preprocess bisulfite sequencing data so that methylation information from one cytosine is "expanded" to neighboring cytosine
#takes in a path to a directory	or a file as arguments
import csv
import sys


# biseq_files = ['meth_seq_CLLE_cds.bed', 'meth_seq_MALY_cds.bed', 'meth_seq_PBCA_cds.bed']
#biseq_files = ['normal_colon_WGBS_processed.bed']
directory = '/data/khandekara2/normal_WGBS/processed_data'
for file in os.listdir(directory):
    if file.endswith(".cds")
        with open(file, 'r') as f, open(file.split('.')[0] + file.split('.')[1] + file.split('.')[2] + file.split('.')[3] + file.split('.')[5] , 'w') as csvout:
            reader = csv.reader(f, delimiter='\t')
            writer = csv.writer(csvout, delimiter='\t')
            rownum = 0
            for row in reader:
                writer.writerow([row[0], row[1], row[2], row[3], row[4]]) #chrom, start, stop, id, methylation ratio
                writer.writerow([row[0], int(row[1]) + 1, int(row[2]) + 1, row[3], row[4]]) #chrom, start, stop, id, methylation ratio


