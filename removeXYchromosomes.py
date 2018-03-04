import csv

with open('normal_colon_WGBS.bed', 'r') as f, open('normal_colon_WGBS_processed.bed', 'w') as csvout:
	reader = csv.reader(f, delimiter='\t')
	writer = csv.writer(csvout, delimiter='\t')
	for row in reader:
		if not row[0] == 'chrX' and not row[0] == 'chrY' :
			writer.writerow(row)