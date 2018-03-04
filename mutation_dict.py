# create a dictionary that maps the genomic location of a CpG mutation to a list of all the samples that that mutation occurred in

import pickle
import pandas as pd
import os

#common_samples = pickle.load(open('common_samples.pickle', 'rb'))

for file in os.listdir("/data/khandekara2/mutation_data"):
        if file.endswith(".noDuplicates.bed.tsv"):
                cancer_type = file.split('_')[0]
                mutation_dict = {} #key is chromosomal location, value is list with id's of all samples that the the mutation occurred in
                df = pd.read_csv(file, sep='\t')
                #add location column
                locations = []
                for chrom, start, stop in zip(list(df['chromosome']), list(df['chromosome_start']), list(df['chromosome_end'])):
                        locations.append((chrom, start, stop))
                df['location'] = locations
                #populate dictionary
                for loc, sample in zip(df.location, df.id):
                        coords = (str(loc[0]), int(loc[1]), int(loc[2]))
                        if loc not in mutation_dict:
                                mutation_dict[coords] = []
                        mutation_dict[coords].append(sample)

                #pickle dictionary
                pickle.dump(mutation_dict, open(cancer_type + '_mutation_dict.pickle', 'wb'))