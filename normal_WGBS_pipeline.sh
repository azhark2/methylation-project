convert2bed -i wig < $f > "$f".bed

FILES=/data/khandekara2/normal_WGBS/processed_data/methylation_mutation/*.non_cds.noSNPS.fixed
EXTENSION=.fixed

for f in $FILES
do

        sed -n "p;N;" $f > temp.temp #delete every other line
        mv temp.temp $f$EXTENSION
        bedtools slop -i $f$EXTENSION -g /data/khandekara2/bed_CpGs/hg19.fa.fai -l 0 -r 1 > temp.temp #adjust coordinates to include both CpG's in dyad
        #mv temp.temp $f$EXTENSION
        bedtools getfasta -fi /data/khandekara2/bed_CpGs/hg19.fa -bed $f -bedOut > temp.temp #make sure that these are in fact CpG's
        mv temp.temp $f
        python check_CpG_context.py #filter out any entries that are not CpG's

done
