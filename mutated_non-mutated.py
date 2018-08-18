#script to split file which contains methylation status of all cytosines in MALY into those that are mutated and non-mutated
import csv
#first we need to create a set that contains tuples of all known mutations
mutations = {}
with open('/data/khandekara2/maly_hotspots/run2/all_overlaps_single.bed', 'r') as csvin:
    reader = csv.reader(csvin, delimiter='\t')
    for row in reader:
        chrom = str(row[0])
        start = int(row[1])
        stop = int(row[2])
        id = str(row[3])
        total_reads = row[9]
        variant_reads = row[10]
        annotation = row[11]
        mutations[(chrom, start, stop, id)] = (total_reads, variant_reads, annotation)


with open ('all_single_2.sorted', 'r') as f, open('all_single.mutated_2.sorted', 'w') as csvout1:
    writer1 = csv.writer(csvout1, delimiter='\t')
    #writer2 = csv.writer(csvout2, delimiter='\t')
    for line in f:
        row = line.strip().split('\t')
        chrom = str(row[0])
        start = int(row[1])
        stop = int(row[2])
        id = str(row[3])
        if (chrom, start, stop, id) in mutations.keys():
            writer1.writerow(row + [mutations[(chrom, start, stop, id)][0], mutations[(chrom, start, stop, id)][1], mutations[(chrom, start, stop, id)][2]])
