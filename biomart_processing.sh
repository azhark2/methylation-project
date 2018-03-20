Steps to process protein coding gene file from Biomart

#rearrange columns so file is in bed format
awk 'BEGIN {FS=OFS="\t"} {print $5,$6,$7,$1,$2,$3,$4,$8,$9,$10,$11,$12,$13}' mart_export.txt > temp.bed

#add chr to first column so file is in bed format
awk 'BEGIN {FS=OFS="\t"} { $1 = "chr" $1; print }' temp.bed.bed > temp.bed

#cut out irrelevant columns
cut -f11 --complement temp.bed > temp2.bed


#use liftOver to convert from GRCH 38 to GRCH 37
