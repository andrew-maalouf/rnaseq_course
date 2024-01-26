#!/usr/bin/env bash

#SBATCH --job-name="gtf2bed"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall

GTF_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/merge_assembly/stringtie_merged.gtf
output_file_1=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/full_bed_file.bed
output_file_2=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/novel_bed_file.bed
output_file_3=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/transcript_bed_file.bed
output_file_4=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/five_prime_bed_file.bed
output_file_5=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/three_prime_bed_file.bed

#module load UHTS/Analysis/bedops/2.4.40

#1: Generate bed file, manually not using gtf2bed, for aesthetic reasons

#this will result in bed file where 1st column is chr number, then start position, then end position, then gene name (if not available add gene_id), then score/1000 (the higher the better), then strandness, then transcript_id
#also semicolons and quotes will be removed
#also, i am just considering lines that start with chr and not GL
awk 'BEGIN{OFS="\t"} $1 ~ /^chr/ && $3=="transcript" {for(i=1;i<=NF;i++) gsub(/[";]/, "", $i); print $1, $4-1, $5, ($14 != "") ? $14 : $10, $6, $7, $12}' "$GTF_DIR" | sort -k1,1 -k2,2n > "$output_file_1"


#2: sort bed file using bedtools sort (mine is sorted)

#3: Generate bed file from step 1 where I only keep novel transcripts aka transcript id $7 starting with MSTRG
grep 'MSTRG' $output_file_1 > $output_file_2


#4: Generate bed file from step 1 where I only keep annotated transcripts aka transcript id $7 starting with ENST
grep 'ENST' $output_file_1 > $output_file_3


#5: Generate bed file from step 1 where I consider positive strand and negative strand to create a window for start position (5' end) -50 nucleotides +50 nucleotides
#then add tab space before columns, otherwise bedtools will throw error
awk ' {if($6 == "+") print $1, $2-50, $2+50, $4, $5, $6, $7; else print $1, $3-50, $3+50, $4, $5, $6, $7}' $output_file_1 | tr ' ' '\t' > $output_file_4



#6: Generate bed file from step 1 where I consider positive and negative strands to create a window for end position (3' end) -50 nucleotides +50 nucleotides
#then add tab space before columns, otherwise bedtools will throw error
awk ' {if($6 == "-") print $1, $2-50, $2+50, $4, $5, $6, $7; else print $1, $3-50, $3+50, $4, $5, $6, $7}' $output_file_1 | tr ' ' '\t' > $output_file_5
