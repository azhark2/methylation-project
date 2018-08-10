#first step is to get the methylation ratio and distance in bp of the closest CpG for every CpG site assayed
#this step is a bit tricky because since this is single base resolution data(at least for most sites) we have to avoid getting the cytosine immediately next to the one in question

#break up large MALY file with all 22 samples into 22 files each containing data for each sample


#sort files by genomic location


#script to get neighboring methylation ratio and distance

#now combine all the 22 sample files into one big file again

#now we need to create a file that contains the average methylation ratio and standard deviation over the 22 samples at each CpG site

#next we use bedtools intersect command to create one file with all the features

#now we separate all the sites based on whether they were mutated or not, the non-mutated sites will be our training data, and the mutated sites will be our out of sample test set

#randomly shuffle the training data and create a validation set

#break up training data into sites that have at least 30 reads and least 60 reads

#now all the input files for training and testing the random forest model are ready





bedtools intersect -a output.bed -b thymus.bed -wa -wb > temp.bed #get normal ratio
cut -f11,12,13 --complement temp.bed > out.bed
bedtools intersect -a out.bed -b /data/khandekara2/cancer_WGBS/raw_data/all_CpGs_MALY_single_averaged.bed.sorted -wa -wb > temp.bed #get tumor average and standard deviation
cut -f8,12,13,14 --complement temp.bed > out.bed
#add header
shuf -n 2000000 out.bed > subset.bed #randomly select 9 million lines to make data smaller and easier to work with

#mutated sites
bedtools intersect -a all_single.mutated_2.sorted -b thymus.bed -wa -wb > temp.bed
cut -f13,14,15 --complement temp.bed > out.bed
bedtools intersect -a out.bed -b /data/khandekara2/cancer_WGBS/raw_data/all_CpGs_MALY_single_averaged.bed.sorted -wa -wb > temp.bed
cut -f14,15,16 --complement temp.bed > out.bed
#chromosome     start   stop    id	methylation_ratio	methylated_reads        unmethylated_reads	neighbor_ratio  distance        total_reads     variant_reads   annotation

#non-mutated sites
