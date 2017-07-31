import csv

#create a bed file with the genomic locations of all bases before a CpG site in in the cds region of the human genome
files = ['meth_seq_CLLE_cds.tsv', 'meth_seq_MALY_cds.tsv', 'meth_seq_PBCA_cds.tsv']
for file in files:
    with open(file, 'r') as f, open('prev_base.bed', 'w') as csvout:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(csvout, delimiter='\t')
        for row in reader:
            chrom = row[0]
            coord = int(row[1]) - 1
            writer.writerow([chrom, coord, coord])






