import csv

with open('temp.fa', 'r') as f, open('meth_seq_MALY_processed.bed', 'w') as csvout:
    reader = csv.reader(f, delimiter='\t')
    writer = csv.writer(csvout, delimiter='\t')
    for row in reader:
        if row[7].upper() == 'G':
            writer.writerow([row[0], int(row[1]) - 1, int(row[2]) - 1, row[3], row[4], row[5], row[6]])
        else:
            writer.writerow([row[0], int(row[1]), int(row[2]), row[3], row[4], row[5], row[6]])


# with open('meth_seq_MALY_shifted.bed', 'r') as f, open('meth_seq_MALY_cds.bed', 'w') as csvout:
#     reader = csv.reader(f, delimiter='\t')
#     writer = csv.writer(csvout, delimiter='\t')
#     for row in reader:
#         writer.writerow([row[0], row[1], row[2], row[3], row[4]])



