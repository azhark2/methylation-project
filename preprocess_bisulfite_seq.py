#preprocesses whole genome bisulfite sequencing files from ICGC by filtering out unneccesary columns and transforming to bed format
import csv
import pickle

#transform bisulfite sequencing files into bed format
biseq_files = ['meth_seq_MALY.tsv', 'meth_seq_PBCA.tsv', 'meth_seq_BOCA.tsv', 'meth_seq_CLLE.tsv']
for file in biseq_files:
    with open(file, 'r') as f, open(file.split('_')[2].split('.')[0] + '.bed', 'w') as csvout:
        reader = csv.reader(f, delimiter='\t')
        csvout = csv.writer(csvout, delimiter='\t')
        rownum = 0
        for row in reader:
            if rownum > 0:
                csvout.writerow(['chr' + str(row[6]), int(row[7]), int(row[8]), row[4], row[11], row[12], row[13]]) #chrom, start, stop, submitted_sample_id, methylation_ratio
            rownum += 1
