sed -i '1i chromosome\tstart\tstop\tid\tmethylation_ratio\ttotal_reads\tmethylated_reads\tfasta' MALY_cds.bed
python bedtools_prep.py
mv *_overlaps.bed maly_hotspots
source combine_overlaps.sh
cut -f11,12,13,14 --complement all_overlaps.bed > temp.bed
mv temp.bed maly_overlaps.bed


#retrieve the gene and other information for "overlaps"
bedtools intersect -a MALY_overlaps.tsv -b mart_export.bed -wa -wb > test.bed
cut -f1,2,3,4,5,6,7,8,9,11,13,15,19,20,21,22  test.bed > test2.bed

#retrieve normal mehtylation ratio for overlaps
