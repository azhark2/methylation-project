#transform Ilumina450k manifest into bed format
import csv
file = 'HumanMethylation450_15017482_v1-2.csv'
with open(file, 'r') as f, open('manifest.bed', 'w') as csvout:
    reader = csv.reader(f)
    writer = csv.writer(csvout, delimiter='\t')
    rownum = 0
    for row in reader:
        if rownum > 7 and not row[11] and not row[12] and not row[0]: # check to make sure no fields are empty
            writer.writerow(['chr' + str(row[11]), row[12], row[12], row[0]]) #chr, start, stop, IlimnID
        rownum += 1

