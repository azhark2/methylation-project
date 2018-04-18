#write out as bed file
python preprocess_bisulfite_seq.py

#adjust coordinates
bedtools shift -i meth_seq_PBCA.bed -g /data/khandekara2/bed_CpGs/hg19.fa.fai -s -1 > meth_seq_PBCA_processed.bed

#expand to include both cytosines in dyad
bedtools slop -i meth_seq_MALY.bed -g /data/khandekara2/bed_CpGs/hg19.fa.fai -l 0 -r 1

#check that adjustment of coordinates worked correctly
bedtools getfasta -fi /data/khandekara2/bed_CpGs/hg19.fa -bed meth_seq_BOCA_processed.bed -bedOut > temp.fasta

#filter for CpG's in protein coding regions
bedtools intersect -a meth_seq_BOCA_processed.bed -b cds.bed > meth_seq_BOCA_cds.bed

bedtools getfasta -fi hg19.fa -bed meth_seq_MALY_processed.bed -bedOut > temp.bed

MALY:

bedtools slop -i meth_seq_MALY.bed -g /data/khandekara2/bed_CpGs/hg19.fa.fai -l 0 -r 2
bedtools shift -i temp.bed -g /data/khandekara2/bed_CpGs/hg19.fa.fai -s -1 > meth_seq_MALY_processed.bed
bedtools getfasta -fi /data/khandekara2/bed_CpGs/hg19.fa -bed meth_seq_MALY_processed.bed -bedOut > temp.bed
python fix_MALY.py
python remove_duplicates.py MALY_cds.bed
mv MALY_cds.bed.noDuplicates MALY_cds.bed


#convert2bed < $f > $f.bed
