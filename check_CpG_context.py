#check to make sure normal WGBS files have ratios only for cytosines in CpG contexts(like the cancer WGBS)
import csv
import os
for file in os.listdir('.'):
    if file.endswith('.wig.bed'):
        with open(file + '_CpGs', 'w') as csvout, open(file, 'r') as f:
            count = 0
            writer = csv.writer(csvout, delimiter='\t')
            reader = csv.reader(f, delimiter='\t')
	        for row in reader:
                if row[7].upper() == 'CG':
                    writer.writerow(row)
            # print (file + ' ' + str(count))
