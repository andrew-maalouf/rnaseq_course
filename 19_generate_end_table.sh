#!/usr/bin/env bash

#SBATCH --job-name="table"
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2
#SBATCH --partition=pall

#list of used files
CODING_POTENTIAL=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/end_tables/required_files/bed_with_coding_potential.bed
FIVE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/end_tables/required_files/five_prime_overlap.bed
THREE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/end_tables/required_files/three_prime_overlap.bed
BED=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/end_tables/required_files/full_bed_file.bed
TRANSCRIPT_LEVEL=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/end_tables/required_files/transcript_level_results.tsv
INTERGENIC=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/end_tables/required_files/intergenic_novel.bed

#Transcript level results has 274000+ lines while Bed file has 273000+ lines
#Keep lines that are common based on transcript_id (transcript_id is $1 in results table and $7 in bed)
awk 'NR==FNR{a[$7]; next} $1 in a' $BED $TRANSCRIPT_LEVEL > reduced_transcript_level_results.tsv

#I am only interested in novel transcripts, hence I will extract from "reduced results" only my novel transcripts (it will be further reduced)
grep "MSTRG" reduced_transcript_level_results.tsv > novel_transcript_level_results.tsv

#delete previous table
rm reduced_transcript_level_results.tsv

#keep in novel transcript results table only these columns in this order: 
#transcript_id ($1); gene_id ($2); biotype ($3); qval ($5); log2fc ($6); mean_obs ($8)
awk -v OFS='\t' '{print $1, $2, $3, $5, $6, $8}' novel_transcript_level_results.tsv > novel_results.tsv
rm novel_transcript_level_results.tsv

#add coding potential for each transcript_id $1 in novel_results.tsv
#transcript_id is $7 in bed with coding potential
#coding potential is $8 in bed with coding potential
awk -v OFS='\t' 'NR==FNR{transcript_id[$7]=$8; next} {print $0, transcript_id[$1]}' $CODING_POTENTIAL novel_results.tsv > merged_file.txt

#rename merged file and delete old table
rm novel_results.tsv
mv merged_file.txt novel_results.tsv

#add column with "X" if this transcript overlaps with 5'
awk 'NR==FNR{transcript_id[$7]; next} {if ($1 in transcript_id) print $0 "\tX"; else print $0 "\t-"}' $FIVE_PRIME novel_results.tsv > merged_file.txt

#add column with "X" if this transcript overlaps with 3'
awk 'NR==FNR{transcript_id[$7]; next} {if ($1 in transcript_id) print $0 "\tX"; else print $0 "\t-"}' $THREE_PRIME merged_file.txt > merged_file_2.txt

#delete old table and intermediary merged file the rename end merged file
rm novel_results.tsv
rm merged_file.txt
mv merged_file_2.txt novel_results.tsv

#add column with "X" if intergenic
awk 'NR==FNR{transcript_id[$7]; next} {if ($1 in transcript_id) print $0 "\tX"; else print $0 "\t-"}' $INTERGENIC novel_results.tsv > merged_file.tsv

#re-ordering columns and adding header
awk 'BEGIN {OFS="\t"; print "transcript_id", "gene_id", "log2fc", "qval", "five_prime", "three_prime", "intergenic", "mean_obs", "coding_potential", "biotype"} {print $1, $2, $5, $4, $8, $9, $10, $6, $7, $3}' merged_file.tsv > end_table.tsv
rm merged_file.tsv
rm novel_results.tsv

#extract transcripts with differential significance
awk -F'\t' 'NR==1 || ($4 < 0.05)' end_table.tsv > diff_signif.tsv

#get ideal transcripts
#criteria: 5' 3' and intergenic + qval < 0.05
awk '$5 == "X" && $6 == "X" && $7 == "X" {print}' diff_signif.txt > ideal.tsv

#extract transcripts with differential significance + not protein coding
awk -F'\t' 'NR==1 || ($9 < 0.364)' diff_signif.tsv > lnc_diff_signif.tsv

#from lnc diff signif, keep transcripts that are intergenic
#these can be considered for candidates if log2fc shows a differential biologic expression
#I will use log2fc + or - 0.99 as threshold just by pure choice
awk -F'\t' 'NR==1 || ($7 == "X")' lnc_diff_signif.tsv > lnc_diff_signif_intergenic.tsv
awk -F'\t' 'NR==1 || ($3 < -0.99 || $3 > 0.99)' lnc_diff_signif_intergenic.tsv > possible_candidates.tsv

#look for candidates based on mean_obs
#sort by mean_obs $8. keep in mind that there is header which is now the last line
awk '{print $0}' lnc_diff_signif.tsv | sort -k8,8rn > sorted_lnc_diff_signif.tsv

