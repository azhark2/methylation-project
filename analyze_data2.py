from __future__ import division
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
import pickle


# create dictionary of cg_id --> gene name, (chromosome, chromosome coordinate), gene region
def create_master_dict(manifest_location):
    # manifest_df = pd.read_csv(manifest_location, header=7)
    # manifest_df = manifest_df.dropna(subset=['UCSC_RefGene_Name'])
    #manifest_df.to_pickle('manifest.pickle')
    manifest_df = pd.read_pickle('manifest.pickle')
    genes = manifest_df.UCSC_RefGene_Name
    chrs = manifest_df.CHR
    chr_locations = manifest_df.MAPINFO
    cg_ids = manifest_df.IlmnID
    gene_regions = manifest_df.UCSC_RefGene_Name
    master_dict = {}
    for cg_id, gene, chr, chr_location, gene_region in zip(cg_ids, genes, chrs, chr_locations, gene_regions):
        master_dict[cg_id] = (gene, (chr, int(chr_location)), gene_region)
    return master_dict

# read in csv data file
def read_csv(filename):
    tp = pd.read_csv(filename, sep='\t', usecols=[3, 7, 8], iterator=True, chunksize=1000, low_memory=False, compression = 'gzip')
    df = pd.concat(tp, ignore_index=True)
    return df


# plot distribution of all beta-values for cancer type
def plot_betaValueDistribution(f, all_beta_values):
    mlist = list(df.iloc[:, 2])
    all_beta_values = all_beta_values + mlist
    plt.hist(mlist, bins=10)
    plt.xticks(np.arange(0, 1.1, 0.1))
    plt.title("Distribution of Beta Values for " + f[len(f) - 14:-7])
    plt.xlabel("Beta Values")
    plt.ylabel("Frequency")
    plt.show()
    plt.savefig(f[len(f) - 14:-6] + 'png')
    plt.clf()





# count total methylated sites in file(beta-value > 0.8)
def count_total_methylated_sites():
    methylated = df.loc[df.iloc[:, 2] >= 0.8]
    me_only = list(methylated.iloc[:, 1])

    return len(me_only)

# count methylated sites in intragenic regions(sites that fall into coding regions of genes)
# create dictionary that maps genes --> # of methylated sites
# create dictionary that maps genes --> # of total sites assayed
def count_gene_methylated_sites():
    cgs = list(df.iloc[:, 1])
    beta_values = list(df.iloc[:, 2])
    gene_methylated = 0
    for cg, beta_value in zip(cgs, beta_values):
        if cg in master_dict:
            gene_tuple = master_dict[cg]
            gene = gene_tuple[0].split(';')
            gene = gene[0]
            if gene not in gene_dict2:
                gene_dict2[gene] = 0
            gene_dict2[gene] += 1
            if beta_value >= 0.8:
                gene_methylated += 1
                if gene not in gene_dict:
                    gene_dict[gene] = 0
                gene_dict[gene] += 1
            if gene not in gene_beta_values:
                gene_beta_values[gene] = []
            gene_beta_values[gene].append(beta_value)

    return gene_methylated


#main

# files = ['https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PAEN-AU/meth_array.PAEN-AU.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BLCA-US/meth_array.BLCA-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-US/meth_array.BRCA-US.tsv.gz',
#        'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/CESC-US/meth_array.CESC-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/CLLE-ES/meth_array.CLLE-ES.tsv.gz',
# 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/GBM-US/meth_array.GBM-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LAML-US/meth_array.LAML-US.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/UCEC-US/meth_array.UCEC-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/STAD-US/meth_array.STAD-US.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-US/meth_array.PRAD-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-CA/meth_array.PRAD-CA.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PBCA-DE/meth_array.PBCA-DE.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/OV-US/meth_array.OV-US.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/OV-AU/meth_array.OV-AU.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/HNSC-US/meth_array.HNSC-US.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/READ-US/meth_array.READ-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/THCA-US/meth_array.THCA-US.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LGG-US/meth_array.LGG-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LIHC-US/meth_array.LIHC-US.tsv.gz',
#          'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/GBM-US/meth_array.GBM-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LAML-US/meth_array.LAML-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/KIRP-US/meth_array.KIRP-US.tsv.gz',
#          'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/KIRC-US/meth_array.KIRC-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/SKCM-US/meth_array.SKCM-US.tsv.gz']
#
# files = ['https://dcc.icgc.org/api/v1/download?fn=/current/Projects/OV-AU/meth_array.OV-AU.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/HNSC-US/meth_array.HNSC-US.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/READ-US/meth_array.READ-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/THCA-US/meth_array.THCA-US.tsv.gz',
#              'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LGG-US/meth_array.LGG-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LIHC-US/meth_array.LIHC-US.tsv.gz',
#          'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/GBM-US/meth_array.GBM-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LAML-US/meth_array.LAML-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/KIRP-US/meth_array.KIRP-US.tsv.gz',
#          'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/KIRC-US/meth_array.KIRC-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/SKCM-US/meth_array.SKCM-US.tsv.gz']

files = ['https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PAEN-AU/meth_array.PAEN-AU.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/GBM-US/meth_array.GBM-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/SKCM-US/meth_array.SKCM-US.tsv.gz'  ]
total_methylated = 0
total_gene_methylated = 0
#create list of Cosmic genes implicated in cancer
#df1 = pd.read_csv('/Users/khandekara2/Documents/methylationProject/01_data/Census_allThu Feb  9 18_47_50 2017.csv', header = 0, usecols =[0]) #Cosmic cancer genes
df1 = pd.read_csv('Census_allThu Feb  9 18_47_50 2017.csv', header = 0, usecols =[0]) #Cosmic cancer genes
cosmicGenes = set(df1.iloc[:, 0]) #use set because they are much faster for membership testing(O(1) as opposed to O(n))

gene_dict = {}  #maps genes to the total number of methylated sites in that gene
gene_dict2 = {} #maps genes to the total number of sites assayed for methylation in that gene
gene_beta_values = {} # maps genes to a list of beta values corresponding to that gene
# manifest_location = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv'
# master_dict = create_master_dict(manifest_location)
# pickle.dump(master_dict, open("master_dict.pickle", 'wb'))
master_dict = pickle.load( open( "master_dict.pickle", "rb" ))
print ('Finished with manifest')
all_beta_values = []
filenames = ['https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-US/meth_array.PRAD-US.tsv.gz', 'https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-CA/meth_array.PRAD-CA.tsv.gz']
siteCount = {} #dictionary that maps cancer type to (# of methylated sites, # of gene methylated sites)
sample_count = []
for f in files:
    try:
        df = pd.read_pickle(f[len(f) - 14:-6] + 'pickle')
        print ('Successfully loaded pickled data frame of' %s(f[len(f) - 14:-6]))
    except:
        df = read_csv(f)
        print("finished reading" + f[len(f) - 14:-6])
    plot_betaValueDistribution(f, all_beta_values)
    numMethylated = count_total_methylated_sites()
    numSamples = df['icgc_sample_id'].nunique()
    sample_count.append(numSamples)
    print ('Number of samples for ' + f[len(f) - 14:-7] + ': %d' % (numSamples))

    print('Number of methylated sites per sample for ' + f[len(f) - 14:-7] + ': %d' % (int (numMethylated / numSamples)))
    total_methylated += numMethylated
    geneMethylated = count_gene_methylated_sites()
    print('Number of methylated sites in coding regions per sample for ' + f[len(f) - 14:-7] + ': %d' % (geneMethylated / numSamples))
    total_gene_methylated += geneMethylated

    #take care of the problem that some cancer types are split into multiple files
    if f[len(f) - 14:-10] in siteCount:
        siteCount[f[len(f) - 14:-10]][0] += numMethylated
        siteCount[f[len(f) - 14:-10]][1] += geneMethylated
    else:
        siteCount[f[len(f) - 14:-10]]= [numMethylated, geneMethylated]


#transform counts into proportions
for key, value in gene_dict.items():
    gene_dict[key] = value / gene_dict2[key]

print ('Number of total methylated sites across all cancer types: %d ' %(total_methylated))
print ('Number of total methylated sites across all cancer types that fall into the coding regions of genes: %d' %(total_gene_methylated))
# print (list(gene_dict.items())[:5])
pickle.dump(gene_dict, open("gene_dict.pickle ", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(gene_dict2, open("gene_dict2.pickle ", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(gene_beta_values, open("gene_beta_values.pickle", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(all_beta_values , open("all_beta_values .pickle", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(siteCount, open("siteCount.pickle", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(sample_count, open('sample_count.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)



#plot distribution of all beta values
plt.hist(all_beta_values, bins=10)
plt.xticks(np.arange(0, 1.1, 0.1))
plt.title("Distribution of ALL Beta Values")
plt.xlabel("Beta Values")
plt.ylabel("Frequency")
plt.savefig('all_beta_values_distribution.png')
plt.show()

# print(list(master_dict.items())[-10:])

topCosmicGenes = []
topCosmicCounts = []
sorted_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse = True) #list of tuples of genes and the proportion of methylated sites(from highest to lowest)
genes = [pair[0] for pair in sorted_genes]
proportion_methylated = [pair[1] for pair in sorted_genes]


for gene, count in zip(genes, proportion_methylated):
    if gene in cosmicGenes:
        topCosmicGenes.append(gene)
        topCosmicCounts.append(count)

#plot top 20 genes with the greatest proportion of methylated sites
plt.title("Top 20 genes with the greatest proportion of methylated sites")
plt.barh(np.arange(20), proportion_methylated[:20], align='center')
plt.yticks(np.arange(20), genes[:20])
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Proportion of Methylated Sites')
plt.ylabel('Genes')
plt.gca().invert_yaxis()
plt.savefig('top_methylated_genes.png')
# plt.show()


#plot top 20 COSMIC CANCER genes with the greatest proportion of methylated sites
plt.title("Top 20 COSMIC CANCER genes with the greatest proportion of methylated sites")
plt.barh(np.arange(20), topCosmicCounts[:20], align='center')
plt.yticks(np.arange(20), topCosmicGenes[:20])
plt.xlabel('Proportion of Methylated Sites')
plt.ylabel('COSMIC Cancer Genes')
plt.gca().invert_yaxis()
# plt.show()
plt.savefig('top_methylated_COSMIC_genes.png')


#plot each cancer type and how many methylated sites there are
cancer_types = []
total_sites = []
gene_sites = []
for cancer_type, count in siteCount.items():
    cancer_types.append(cancer_type)
    total_sites.append(count[0])
    gene_sites.append(count[1])

plt.title("# of Methylated Sites per Cancer Type")
plt.barh(np.arange(len(cancer_types)), total_sites, align='center')
plt.yticks(np.arange(len(cancer_types)), cancer_types)
plt.savefig('cancer_type_methylated_sites.png')

#plot each cancer type and how many methylated sites in protein coding regions there are
plt.title("# of Gene Methylated Sites per Cancer Type")
plt.bar(np.arange(len(cancer_types)), gene_sites, align='center')
plt.xticks(np.arange(len(cancer_types)), cancer_types, rotation=90)
plt.savefig('cancer_type_gene_methylated_sites.png')
# plt.show()
