#shell script to produce regressions plot

#first step is to process the WGBS data from normal tissues
FILES=/panfs/pan1.be-md.ncbi.nlm.nih.gov/methylation-mutations/normal_WGBS/raw_data/*.wig #raw data files in wig format
EXTENSION=.processed
for f in $FILES
do

	convert2bed -i wig < $f > temp.temp #convert wig files to bed files
  mv temp.temp $f$EXTENSION
  sed -n "p;N;" $f$EXTENSION > temp.temp #delete every other line
  mv temp.temp $f$EXTENSION
  bedtools slop -i $f$EXTENSION -g /data/khandekara2/bed_CpGs/hg19.fa.fai -l 0 -r 1 > temp.temp #adjust coordinates to include both CpG's in dyad
  mv temp.temp $f$EXTENSION
  bedtools getfasta -fi /data/khandekara2/bed_CpGs/hg19.fa -bed $f$EXTENSION -bedOut > temp.temp #make sure that these are in fact CpG's by getting fasta
  mv temp.temp $f$EXTENSION
  python check_CpG_context.py #filter out any entries that are not CpG's
  bedtools intersect -a $f$EXTENSION -b /data/khandekara2/bed_CpGs/single_SNPS.bed -v > temp.temp #filter out all SNPS
  mv temp.temp $f$EXTENSION

done

#second step is to process the mutation data from ICGC
python preprocess_mutations.py ##filters all files for only C>T or G>A mutations and removes noDuplicates
FILES=/data/khandekara2/mutation_data/raw_data/*.noDuplicates
for f in $FILES
do
  sed -i '1i #chromosome\tchromosome_start\tchromosome_end\tid\tmutated_from_allele\tmutated_to_allele\ttotal_read_count\tmutant_allele_read_count\tconsequence_type\tgene_affected' $f #add header in place
  #adjust coordinates
  bedtools slop -i $f -g hg19.fa.fai -l 1 -r 0 > temp.bed #adjust coordinates
  mv temp.bed $f$BED
done
