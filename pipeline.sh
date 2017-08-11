#write out as bed file
python preprocess_bisulfite_seq.py 


#adjust coordinates
bedtools shift -i meth_seq_PBCA.bed -g hg19.fa.fai -s -1 > meth_seq_PBCA_processed.bed 

#check that adjustment of coordinates worked correctly
bedtools getfasta -fi hg19.fa -bed meth_seq_BOCA_processed.bed -bedOut > temp.fasta

#filter for CpG's in protein coding regions
bedtools intersect -a meth_seq_BOCA_processed.bed -b cds.bed > meth_seq_BOCA_cds.bed


bedtools getfasta -fi hg19.fa -bed meth_seq_MALY_processed.bed -bedOut > temp.bed

MALY:

bedtools slop -i meth_seq_MALY.bed -g hg19.fa.fai -l 0 -r 1
bedtools shift -i temp.bed -g hg19.fa.fai -s -1 > meth_seq_MALY_processed.bed 