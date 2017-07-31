# Goal: Figure out if there are overlapping samples between bisulfite sequencing data and somatic mutation data

import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np
import pickle


wgbs_files = ['meth_seq_MALY_cds.tsv', 'meth_seq_CLLE_cds.tsv', 'meth_seq_PBCA_cds.tsv']
ssm_files = ['download?fn=%2Fcurrent%2FProjects%2FMALY-DE%2Fsimple_somatic_mutation.open.MALY-DE.tsv', 'download?fn=%2Frelease_23%2FProjects%2FCLLE-ES%2Fsimple_somatic_mutation.open.CLLE-ES.tsv',
             'download?fn=%2Frelease_23%2FProjects%2FPBCA-DE%2Fsimple_somatic_mutation.open.PBCA-DE.tsv']



list_of_coding_cpgs = pickle.load(open('list_of_coding_cpgs.pickle', 'rb'))

cancer_cpgs = set([])
file = 'cancer_cds_cpg.bed'
with open(file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        cancer_cpgs.add((row[0], (row[1]))) #tuple of strings
        cancer_cpgs.add((row[0], (row[2])))
avgBeta = {} #key: cancer type value: dictionary that maps location to avg methylation ratio observed at that location
cancer_overlaps = {} #key: cancer type value: list of lists(overlaps in cds, overlaps in cancer cds)

with open("results.csv", 'w') as csvout:
    actual_ratios = []
    writer = csv.writer(csvout)

    #header
    writer.writerow(['cancer_type', 'sample_id', 'num_sites_covered', 'num_cds_sites_methylated', 'num_cancer_cds_sites_methylated', 'num_non_cancer_cds_sites_methylated', 'num_cds_cpg_sites_mutated', 'num_cancer_cds_cpg_sites_mutated', 'num_non_cancer_cds_cpg_sites_mutated', 'num_overlaps',
                     'num_cancer_overlaps', 'num_non_cancer_overlaps'])


    for wgbs_file, ssm_file in zip(wgbs_files, ssm_files):
        cancer_type = wgbs_file.split('_')[2]
        all_cds_overlaps = [] #all overlaps for particular cancer type
        all_cancer_cds_overlaps = []
        avgBeta[cancer_type] = {}
        biseq = pd.read_csv(wgbs_file, sep='\t')
        mutation = pd.read_csv(ssm_file, sep='\t')
        locations = []
        locations2 = []

        #append extra column to dataframes that specifies genomic location
        for chrom, coord in zip(list(biseq['chromosome']), list(biseq['chromosome_start'])):
            locations.append(('chr' + str(chrom), str(coord)))
        biseq['location'] = locations

        for chrom, coord in zip(list(mutation['chromosome']), list(mutation['chromosome_start'])):
            locations2.append(('chr' + str(chrom), str(coord)))
        mutation['location'] = locations2



        #create dictionary of site > avg methylation ratio
        groupby = biseq['methylation_ratio'].groupby(biseq['location'])
        df = groupby.agg([np.mean])
        df.reset_index(inplace=True)
        for location, mean in zip(df['location'], df['mean']):
            avgBeta[cancer_type][location] = mean

        #get common samples
        biseq_samples = set(biseq['submitted_sample_id'].unique())
        mutation_samples = set(mutation['submitted_sample_id'].unique())
        common_samples = biseq_samples & mutation_samples
        print(len(biseq_samples))
        print(len(mutation_samples))
        print(len(common_samples))

    # Goal: to do a sample by sample comparison of mutated and methylated sites CpG sites and write them out to csv file


        for sample in list(common_samples):
            # subset dataframes for specific samples
            biseq['submitted_sample_id'] = biseq['submitted_sample_id'].str.strip()
            mutation['submitted_sample_id'] = mutation['submitted_sample_id'].str.strip()

            sub_biseq = biseq[biseq['submitted_sample_id'] == sample]
            sub_mutation = mutation[mutation['submitted_sample_id'] == sample]

            methylated = sub_biseq[sub_biseq['methylation_ratio'] > 0.7]

            methylated_cds = set([])  # all methylated sites in cds
            methylated_cancer_cds_cpg = set([])  # all methylated sites in cancer cds
            mutation_cds = set([])
            mutation_cancer_cds = set([])  # all mutated sites in cancer cds cpgs

            siteToRatio = {} #dictionary that maps each site in sample to its methylation ratio
            for location, ratio in zip(sub_biseq['location'], sub_biseq['methylation_ratio']):
                siteToRatio[location] = ratio

            # find methylated sites in
            for chrom, coord in zip(list(methylated['chromosome']), list(methylated['chromosome_start'])):
                location = ('chr' + str(chrom), str(coord))
                methylated_cds.add(location)
                if location in cancer_cpgs:
                    methylated_cancer_cds_cpg.add(location)

            sub_mutation = sub_mutation[
                (sub_mutation.reference_genome_allele == 'C') & (sub_mutation.mutated_to_allele == 'T') | (
                    sub_mutation.reference_genome_allele == 'G') & (
                    sub_mutation.mutated_to_allele == 'A')]
            sub_mutation = sub_mutation[sub_mutation.consequence_type == 'missense_variant']

            asm = sub_biseq[sub_biseq['methylation_ratio'] > 0.4]  # allele specific methylation
            asm_cds = set([])  # all partially methylated sites
            asm_cancer_cds = set([])

            # find partially methylated sites in cds and cancer cds
            for location in sub_biseq['location']:
                asm_cds.add(location)
                if location in cancer_cpgs:
                    asm_cancer_cds.add(location)

            # find mutated CpG sites in cds and cancer cds
            for location in sub_mutation['location']:
                if location in list_of_coding_cpgs:
                    mutation_cds.add(location)
                if location in cancer_cpgs:
                    mutation_cancer_cds.add(location)

            # lists of overlaps
            cds_overlap = set.intersection(asm_cds, mutation_cds)  #mutated and methylated sites in cds


            for overlap in cds_overlap:
                actual_ratios.append(siteToRatio[overlap])

            cancer_cds_overlap = set.intersection(asm_cancer_cds, mutation_cancer_cds) #mutated and methylated sites in cancer cds


            #keep track of all the overlaps in cds of all genes and cds of cancer genes
            all_cds_overlaps.extend(cds_overlap)
            all_cancer_cds_overlaps.extend(cancer_cds_overlap)



            #get beta values for overlaps


            # write results to csv file:
            # num sites covered, num cds sites methylated, num cancer cds sites methylated, num non cancer cds sites methylated,
            # num cds cpg sites mutated, num cancer cds cpg sites mutated, num non cancer cds cpg sites mutated, num overlaps, num cancer overlaps, num non cancer overlaps
            writer.writerow(
                [wgbs_file.split('_')[2], sample, sub_biseq.shape[0], len(methylated_cds), len(methylated_cancer_cds_cpg), len(methylated_cds) -
                 len(methylated_cancer_cds_cpg), len(mutation_cds), len(mutation_cancer_cds),
                 sub_mutation.shape[0] - len(mutation_cancer_cds), len(cds_overlap), len(cancer_cds_overlap),
                 len(cds_overlap) - len(cancer_cds_overlap)])

        cancer_overlaps[cancer_type] = []
        cancer_overlaps[cancer_type].append(all_cds_overlaps)
        cancer_overlaps[cancer_type].append(all_cancer_cds_overlaps)


# pickle.dump(avgBeta, open('WGBS_avgBeta.pickle', 'wb'))
pickle.dump(cancer_overlaps, open('cancer_overlaps_no_threshold.pickle', 'wb'))
pickle.dump(actual_ratios, open('actual_ratios_no_threshold.pickle', 'wb'))




