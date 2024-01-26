#!/usr/bin/env bash

#SBATCH --job-name="stat"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall


MY_FIVE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/five_prime_bed_file.bed
MY_THREE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/three_prime_bed_file.bed
NOVEL_BED=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/novel_bed_file.bed
TRANSCRIPT_BED=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/transcript_bed_file.bed
OUTPUT_INTERGENIC=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/intergenic_novel.bed
FIVE_PRIME_OVERLAP=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/five_prime_overlap.bed
THREE_PRIME_OVERLAP=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/three_prime_overlap.bed
CODING_POTENTIAL=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/coding_potential.bed
WRITTEN_ANSWERS=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/written_answers.txt
IG_NC=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/intergenic_novel_non_coding.bed
IG_C=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/intergenic_novel_coding.bed
N_C=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/novel_coding.bed
REDUCED_FIVE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/just_novel_transcripts_5.bed
REDUCED_THREE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/just_novel_transcripts_3.bed
REDUCED_FIVE_PRIME_OVERLAP=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/just_novel_transcripts_overlap_5.bed
REDUCED_THREE_PRIME_OVERLAP=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/just_novel_transcripts_overlap_3.bed



#QUESTION 3: Intergenic novels
echo "Question 3: How many novel “intergenic” (not overlapping annotated protein coding gene spans) genes have you identified?" > $WRITTEN_ANSWERS
total_novel=$(wc -l $NOVEL_BED | awk '{print $1}' )
intergenic_novel=$(wc -l $OUTPUT_INTERGENIC | awk '{print $1}' )
percentage_intergenic_novel=$(echo " scale=4; ($intergenic_novel / $total_novel) * 100" | bc)
echo "Out of $total_novel considered novels, $intergenic_novel are intergenic, meaning $percentage_intergenic_novel %" >> $WRITTEN_ANSWERS 
echo "" >> $WRITTEN_ANSWERS
#QUESTION 3/2: Intergenic novels and coding potential:
echo "Question 2: What percent of your novel transcripts are protein coding?" >> $WRITTEN_ANSWERS
echo "If coding potential threshold is 0.364, then:" >> $WRITTEN_ANSWERS
intergenic_novel_non_coding=$(wc -l $IG_NC | awk '{print $1}')
intergenic_novel_coding=$(wc -l $IG_C | awk '{print $1}')
percentage_nc=$(echo " scale=4; ($intergenic_novel_non_coding / $intergenic_novel) * 100" | bc)
percentage_c=$(echo " scale=4; ($intergenic_novel_coding / $intergenic_novel) * 100" | bc)
echo "Out of $intergenic_novel intergenic novels, $intergenic_novel_non_coding are non coding, meaning $percentage_nc %" >> $WRITTEN_ANSWERS
echo "and $intergenic_novel_coding are coding, meaning $percentage_c %" >> $WRITTEN_ANSWERS

novel_coding=$(wc -l $N_C | awk '{print $1}')
percentage_coding_novel=$(echo " scale=4; ($novel_coding / $total_novel) * 100" | bc)
echo "Out of $total_novel considered novels, $novel_coding are coding, meaning $percentage_coding_novel %" >> $WRITTEN_ANSWERS 

#Question 1:
echo "" >> $WRITTEN_ANSWERS
echo "Question 1: How good are the 5’ and 3’ annotations of your transcripts?" >> $WRITTEN_ANSWERS
my_total_five=$(wc -l $MY_FIVE_PRIME | awk '{print $1}')
my_total_three=$(wc -l $MY_THREE_PRIME | awk '{print $1}')
overlap_five=$(wc -l $FIVE_PRIME_OVERLAP| awk '{print $1}')
overlap_three=$(wc -l $THREE_PRIME_OVERLAP | awk '{print $1}')
percentage_five=$(echo " scale=4; ($overlap_five / $my_total_five) * 100" | bc)
percentage_three=$(echo " scale=4; ($overlap_three / $my_total_three) * 100" | bc)
echo "$percentage_five % of my 5' annotations overlap" >> $WRITTEN_ANSWERS
echo "$percentage_three % of my 3' annotations overlap" >> $WRITTEN_ANSWERS
#if considering novel transcripts only
my_total_five=$(wc -l $REDUCED_FIVE_PRIME | awk '{print $1}')
my_total_three=$(wc -l $REDUCED_THREE_PRIME | awk '{print $1}')
overlap_five=$(wc -l $REDUCED_FIVE_PRIME_OVERLAP| awk '{print $1}')
overlap_three=$(wc -l $REDUCED_THREE_PRIME_OVERLAP | awk '{print $1}')
percentage_five=$(echo " scale=4; ($overlap_five / $my_total_five) * 100" | bc)
percentage_three=$(echo " scale=4; ($overlap_three / $my_total_three) * 100" | bc)
echo "If considering novel transcripts only:" >> $WRITTEN_ANSWERS
echo "$percentage_five % of my 5' annotations overlap" >> $WRITTEN_ANSWERS
echo "$percentage_three % of my 3' annotations overlap" >> $WRITTEN_ANSWERS