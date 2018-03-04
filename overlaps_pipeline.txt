sed -i '1i chromosome\tstart\tstop\tid\tmethylation_ratio\total_reads\tmethylated_reads' meth_seq_MALY_processed_non_cds.bed
python bedtools_prep.py
mv *.overlaps.bed hotspots3
source combine_overlaps.sh
python remove_duplicates.py all_overlaps_2.bed
cut -f11,12,13,14 --complement all_overlaps_2.bed.noDuplicates > temp.bed
mv temp.bed all_overlaps_2.bed