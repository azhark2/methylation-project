#script to break MALY WGBS file into several files each with a samples information(26 total samples)
#this task is well suited to multi-processing
import csv
from multiprocessing import Pool
import pybedtools
import pandas as pd

def process(sample):
    with open('MALY.fasta', 'r') as csvin, open('MALY_' + sample + '_WGBS.bed', 'w') as csvout:
        reader = csv.reader(csvin, delimiter='\t')
        writer = csv.writer(csvout, delimiter='\t')
        for row in reader:
            if sample == row[3]:
                writer.writerow(row)

if __name__ == '__main__':
    common_samples = ['tumor_4105105', 'tumor_4112512', 'tumor_4119027', 'tumor_4121361', 'tumor_4125240', 'tumor_4133511', 'tumor_4134005', 'tumor_4142267', 'tumor_4158726', 'tumor_4159170', 'tumor_4175837', 'tumor_4177376', 'tumor_4177434', 'tumor_4177856', 'tumor_4182393', 'tumor_4188900', 'tumor_4189200', 'tumor_4189998', 'tumor_4190495', 'tumor_4193278', 'tumor_4194218', 'tumor_4194891'] #samples in MALY that have methylation and mutation data
    pbca_samples = [ f.split('.')[0][4:] for f in os.listdir('/data/khandekara2/bed_CpGs/hotspots3/') if f.startswith('MALY')]
    print (len(common_samples))
    pool = Pool()  # Create a multiprocessing Pool
    pool.map(process, common_samples)
    pool.close()
    pool.join()
# now intersect wgbs files with mutation files to get overlap sites

import csv
from multiprocessing import Pool
import pybedtools
import pandas as pd

def intersect(sample):
    df = pd.read_csv('/data/khandekara2/mutation_data/processed_data/MALY_mutation.bed.all.noDuplicates.tsv', sep='\t')
    sub_mutation = df[df['id'] == sample]
    sub_mutation.to_csv('/data/khandekara2/maly_hotspots/whole_genome/' + sample + '.mutation', sep='\t', index=False, header=False)
    a = pybedtools.BedTool('/data/khandekara2/maly_hotspots/whole_genome/' + sample + '.mutation')
    a.intersect('/data/khandekara2/cancer_WGBS/raw_data/'+ 'MALY_' + sample + '_WGBS.bed', wa=True, wb=True, v=True).saveas('/data/khandekara2/maly_hotspots/whole_genome/' + sample + '_non_mutated.bed')

if __name__ == '__main__':
    common_samples = ['tumor_4105105', 'tumor_4112512', 'tumor_4119027', 'tumor_4121361', 'tumor_4125240', 'tumor_4133511', 'tumor_4134005', 'tumor_4142267', 'tumor_4158726', 'tumor_4159170', 'tumor_4175837', 'tumor_4177376', 'tumor_4177434', 'tumor_4177856', 'tumor_4182393', 'tumor_4188900', 'tumor_4189200', 'tumor_4189998', 'tumor_4190495', 'tumor_4193278', 'tumor_4194218', 'tumor_4194891'] #samples in MALY that have methylation and mutation data
    pool = Pool()  # Create a multiprocessing Pool
    pool.map(intersect, common_samples)
    pool.close()
    pool.join()
