#!/bin/bash

#chromosome     chromosome_start        chromosome_end  submitted_sample_id     mutated_from_allele     mutated_to_allele	consequence_type

python preprocess_mutations.py #includes function for removing duplicates

FILES=/data/khandekara2/mutation_data/raw_data/*.noDuplicates
echo $FILES
TSV=.tsv
BED=.bed

for f in $FILES
do
  	echo $f

        #adjust coordinates
        bedtools slop -i $f -g hg19.fa.fai -l 1 -r 0 > temp.bed #adjust coordinates
        mv temp.bed $f$BED

        #add header back in
        sed '1i chromosome\tchromosome_start\tchromosome_end\tid\tmutated_from_allele\tmutated_to_allele\ttotal_read_count\tmutant_allele_read_count\tconsequence_type\tgene_affected' $f$BED > temp.tsv
        mv temp.tsv $f$TSV

        #bedtools getfasta -fi hg19.fa -bed $f$BED -bedOut > temp.bed
        #mv temp.bed $f$BED

        #sort -k1,1 -k2,2n -k3,3n -u $f > $f.noDuplicates #remove duplicates

done

python /data/khandekara2/mutation_data/mutation_dict.py #use file to create dictionary that maps chromosomal location of mutation to all the samples that mutation occurred in
