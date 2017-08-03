import pybedtools
import pickle
import csv
import os

cancerGeneToCpGDensity = {} #dictionary that maps each Cosmic census cancer gene to its CpG density(# of CpG's in coding region/length of coding region)

a = pybedtools.BedTool('cds.bed')
b = pybedtools.BedTool('all_cds_cpgs_shifted_2.bed.noDuplicates')

#returns total length of cds
def get_cds_length(file):
    if os.stat(file).st_size == 0:
        return 0
    total_length = 0
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            length = int(row[2]) - int(row[1])
            total_length += length
        return total_length


def file_len(fname):
    with open(fname) as f:
        return sum(1 for _ in f)


with open('cancer_genes_hg19.bed', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        with open('temp_file.bed', 'w') as csvout:
            #get all cds corresponding to gene
            writer = csv.writer(csvout, delimiter='\t')
            writer.writerow([row[0], row[1], row[2], row[3]])
        a.intersect('temp_file.bed').saveas('temp2.bed')
        cds_length = get_cds_length('temp2.bed')
        b.intersect('temp2.bed').saveas('temp3.bed')
        num_cpgs = file_len('temp3.bed')
        if cds_length == 0:
            print(row[0], row[1], row[2], row[3])
        else:
            density = num_cpgs / cds_length
            cancerGeneToCpGDensity[row[3]] = density

pickle.dump(cancerGeneToCpGDensity, open('cancerGeneToCpGDensity.pickle', 'wb'))














