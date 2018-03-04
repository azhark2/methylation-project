import csv
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# maly = pd.read_csv('meth_seq_MALY.bed', sep='\t')
# pbca = pd.read_csv('meth_seq_PBCA.bed', sep='\t')
maly_cds = pd.read_csv('meth_seq_MALY_cds.tsv.noDuplicates', sep='\t')
pbca_cds = pd.read_csv('meth_seq_PBCA_cds.tsv.noDuplicates', sep='\t')		


pbca_totalReads = []
maly_totalReads = []
pbca_cds_totalReads = []
maly_cds_totalReads = []

#add a total reads column to each dataframe
# for m, u in zip(pbca['methylated_reads'], pbca['unmethylated_reads']):
#     pbca_totalReads.append(int(m + u))
    
# for m, u in zip(maly['methylated_reads'], maly['unmethylated_reads']):
#     maly_totalReads.append(int(m + u))

for m, u in zip(pbca_cds['methylated_reads'], pbca_cds['unmethylated_reads']):
    pbca_cds_totalReads.append(int(m + u))
    
for m, u in zip(maly_cds['methylated_reads'], maly_cds['unmethylated_reads']):
    maly_cds_totalReads.append(int(m + u))

# pbca['total_reads'] = pbca_totalReads
# maly['total_reads'] = maly_totalReads
pbca_cds['total_reads'] = pbca_cds_totalReads
maly_cds['total_reads'] = maly_cds_totalReads


#plot ratios before applying any threshold
# plt.hist(pbca['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("PBCA Whole Genome, No Read Threshold")
# plt.savefig('PBCA_no_read_threshold.png')
# plt.clf()
#
print len(pbca_cds['methylation_ratio']) 
# plt.hist(pbca_cds['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("PBCA Coding Region, No Read Threshold")
# plt.savefig('PBCA_cds_no_read_threshold.png')
# plt.clf()

# plt.hist(maly['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("MALY Whole Genome, No Read Threshold")
# plt.savefig('MALY_no_read_threshold.png')
# plt.clf()

print len(maly_cds['methylation_ratio']) 
# plt.hist(maly_cds['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("MALY Coding Region, No Read Threshold")
# plt.savefig('MALY_cds_no_read_threshold.png')
# plt.clf()

#plot ratios after applying read threshold of ten
# pbca = pbca[pbca['total_reads'] > 10]
# maly = maly[maly['total_reads'] > 10]
pbca_cds = pbca_cds[pbca_cds['total_reads'] > 10]
maly_cds = maly_cds[maly_cds['total_reads'] >= 10]

# plt.hist(pbca['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("PBCA Whole Genome, Read Threshold = 10")
# plt.savefig('PBCA_read_threshold_ten.png')
# plt.clf()

print len(pbca_cds['methylation_ratio']) 
# plt.hist(pbca_cds['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("PBCA Coding Region, Read Threshold = 10")
# plt.savefig('PBCA_cds_read_threshold_ten.png')
# plt.clf()

# plt.hist(maly['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("MALY Whole Genome, Read Threshold = 10")
# plt.savefig('MALY_read_threshold_ten.png')
# plt.clf()

print len(maly_cds['methylation_ratio']) 
# plt.hist(maly_cds['methylation_ratio'])
# plt.xlabel('methylation ratio')
# plt.ylabel('Frequency')
# plt.title("MALY Coding Region, Read Threshold = 10")
# plt.savefig('MALY_cds_read_threshold_ten.png')
# plt.clf()







