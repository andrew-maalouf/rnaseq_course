#!/usr/bin/env bash

#SBATCH --job-name="analysis"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall

#load modules
module load  UHTS/Analysis/BEDTools/2.29.2
module add SequenceAnalysis/GenePrediction/cpat/1.2.4


FIVE_PRIME_REF=/data/courses/rnaseq_course/lncRNAs/Project1/references/refTSS_v4.1_human_coordinate.hg38.bed
THREE_PRIME_REF=/data/courses/rnaseq_course/lncRNAs/Project1/references/atlas.clusters.2.0.GRCh38.96.bed
MY_FIVE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/five_prime_bed_file.bed
MY_THREE_PRIME=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/three_prime_bed_file.bed
NOVEL_BED=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/novel_bed_file.bed
TRANSCRIPT_BED=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/transcript_bed_file.bed
OUTPUT_INTERGENIC=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/intergenic_novel.bed
FIVE_PRIME_OVERLAP=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/five_prime_overlap.bed
THREE_PRIME_OVERLAP=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/three_prime_overlap.bed
CODING_POTENTIAL=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/coding_potential
WRITTEN_ANSWERS=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/written_answers.txt
BED_FILE=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/full_bed_file.bed
REF_GEN=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/genome/GRCh38.p14.genome.fa
HUMAN_HEXAMER=/data/courses/rnaseq_course/lncRNAs/Project1/references/Human_Hexamer.tsv
LOGIT=/data/courses/rnaseq_course/lncRNAs/Project1/references/Human_logitModel.RData
BED_TO_FASTA=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/bed_to_fasta.fa

#Question 3: finding overlapping regions
#use bedtools intersect with '-v' to keep regions that don't overlap
#these regions are called intergenic
#comparing file a to file b
bedtools intersect -v -a $NOVEL_BED -b $TRANSCRIPT_BED > $OUTPUT_INTERGENIC



#Question 1: finding overlapping 5' transcription starts sites
#window considered in bed file is start position -50 --> start position + 50
# -wa: Write the original entry in A for each overlap
# -s : consider positive strandedness!!!
bedtools intersect -wa -s -a $MY_FIVE_PRIME -b $FIVE_PRIME_REF > $FIVE_PRIME_OVERLAP


#QUESTION 1: finding overlapping 3' sites
#window considered in bed file is start position -50 --> start position + 50
# -wa: Write the original entry in A for each overlap
# -s : consider positive strandedness!!!
bedtools intersect -wa -s -a $MY_THREE_PRIME -b $THREE_PRIME_REF > $THREE_PRIME_OVERLAP


#QUESTION 2: Protein coding potential
# -g: gene file, bed file or fasta file: bed file should be 12 columns -> use fasta -> generate fasta from bed with bedtools
# -o: output file
# -d: logit model
# -x: hexamer table frequency
# if gene file is bed file -> specify reference genome sequence: -r : reference genome sequence

#1st generate fasta file:
# -fo output file name
#  -s force strandedness
# -fi input file name
# -bed bed file


bedtools getfasta -fi $REF_GEN\
 -bed $BED_FILE\
 -s\
 -fo $BED_TO_FASTA

cpat.py -g $BED_TO_FASTA\
 -o $CODING_POTENTIAL\
 -d $LOGIT\
 -x $HUMAN_HEXAMER\
 
