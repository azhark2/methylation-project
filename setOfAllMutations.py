import csv
import pickle
all_mutations = set([]) # set containing tuple of genomic location(chrom, start, stop)
with open('all_mutations.bed', 'r') as csvin:
        reader = csv.reader(csvin, delimiter='\t')
        for row in reader:
                all_mutations.append((row[1], row[2], row[3]))

pickle.dump(all_mutations, open('all_mutations.pickle', 'wb'))


	
