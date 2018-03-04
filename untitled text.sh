#!/bin/bash

#chromosome     chromosome_start        chromosome_end  submitted_sample_id     mutated_from_allele     mutated_to_allele	consequence_type

awk '$1="chr"$1' $1 > temp.bed #add 'chr' to every line in column 1
mv temp.bed $1
tail -n +2 $1 > temp.bed
mv temp.bed $1
tr ' ' '\t' < $1 > temp.bed
mv temp.bed $1
bedtools shift -i $1 -g hg19.fa.fai -s -1 > temp.bed
mv temp.bed $1
bedtools slop -i $1 -g hg19.fa.fai -l 1 -r 0 > temp.bed
mv temp.bed $1
bedtools getfasta -fi hg19.fa -bed $1



awk -F"\t" '!seen[$1, $2, $3, $4]++' $1 > temp.bed
mv temp.bed $1





