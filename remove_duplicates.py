import csv
import sys
#remove duplicates

file = sys.argv[1]

duplicates = set([])
count = 0

with open(file, 'r') as f, open(file + '.noDuplicates', 'w') as csvout:
    reader = csv.reader(f, delimiter='\t')
    writer = csv.writer(csvout, delimiter='\t')
    for row in reader:
        info = (row[0], row[1], row[2], row[3])
        if info not in duplicates:
            duplicates.add(info)
            writer.writerow(row)
        else:
            count += 1