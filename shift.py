import csv

#obtain the coordinates of the neighboring cytosine in CpG dyad
with open('maly_overlaps.bed', 'r') as f, open('maly_overlap_neighbors.bed', 'w') as csvout:
    reader = csv.reader(f, delimiter='\t')
    writer = csv.writer(csvout, delimiter='\t')
    for row in reader:
        if row[13].upper() == 'C':
            writer.writerow([row[0], int(row[1])+1, int(row[2])+1, row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13]])
        else:
            writer.writerow([row[0], int(row[1])-1, int(row[2])-1, row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13]])
