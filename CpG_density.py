from __future__ import division
from multiprocessing import Pool
import pybedtools
import pickle
import csv
import os

cancerGeneToCpGDensity = {} #dictionary that maps each Cosmic census cancer gene to its CpG density(# of CpG's in coding region/length of coding region)
noncancerGeneToCpGDensity = {}

a = pybedtools.BedTool('cds_canonical.bed')
b = pybedtools.BedTool('all_CpGs.bed')

def process(row):
        try:
            	with open(row[3] + '.bed', 'w') as csvout:
                        # get all cds corresponding to gene
                        writer = csv.writer(csvout, delimiter='\t')
                        writer.writerow([row[0], row[1], row[2], row[3]])

                a.intersect(row[3] + '.bed').saveas(row[3] + '_cds.bed')
                file = row[3] + '_cds.bed'

                # get cds length
                cds_length = 0
                if os.stat(file).st_size != 0:
                        with open(file, 'r') as f2:
                                reader2 = csv.reader(f2, delimiter='\t')
                                for row2 in reader2:
                                        length = int(row2[2]) - int(row2[1])
                                        cds_length += length
                print(cds_length)

                # now get number of CpGs in the cds
                b.intersect(row[3] + '_cds.bed').saveas(row[3] + '_cds_CpGs.bed')
                with open(row[3] + '_cds_CpGs.bed') as f3:
                        num_cpgs = sum(1 for _ in f3)

                if cds_length == 0:
                        print(row[0], row[1], row[2], row[3] + ' had cds length of 0')
                else:
                     	if num_cpgs == 0:
                                print(row[0], row[1], row[2], row[3] + ' had 0 CpGs')
                        density = num_cpgs / cds_length
                        print(num_cpgs)
                        noncancerGeneToCpGDensity[row[3]] = density
        except IOError:
                print ("oops file did not exist")

if __name__ == '__main__':
        with open('non_cancer_genes.bed', 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                pool = Pool()  # Create a multiprocessing Pool
                pool.map(process, reader) # process data_inputs iterable with pool
                pickle.dump(noncancerGeneToCpGDensity, open('noncancerGeneToCpGDensity.pickle', 'wb'))
