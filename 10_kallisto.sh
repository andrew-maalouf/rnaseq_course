#!/usr/bin/env bash

#SBATCH --mail-user=andrew.maalouf@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="kallisto"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_error_%j.e

GENOME_FA=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/genome/GRCh38.p14.genome.fa
WORKING_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_expression
ALL_GTF=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/merge_assembly/stringtie_merged.gtf
SAMPLE_FA=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/fastqc/PE

#load modules
module load UHTS/Analysis/kallisto/0.46.0
module add UHTS/Assembler/cufflinks/2.2.1

#create directory and enter it
mkdir /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_expression
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_expression
#create directories for each replicate quantification step
mkdir $WORKING_DIR/H_1_1
mkdir $WORKING_DIR/H_1_2
mkdir $WORKING_DIR/H_1_5
mkdir $WORKING_DIR/P1
mkdir $WORKING_DIR/P2
mkdir $WORKING_DIR/P3

#create fasta file for transcriptome 
gffread -w $WORKING_DIR/transcriptome.fa \
 -g $GENOME_FA \
 $ALL_GTF


#build kallisto index
#--index: name of the file to be generated
# followed by transcriptome reference fasta file to index
kallisto index -i kallisto_genome_index $WORKING_DIR/transcriptome.fa

#kallisto quantification for each paired-end reads
#-i or --index= to specify the index to be used
#-b bootstrap
#-o or --output-dir= directory to write output to
#followed by specifying fasta files of paired-end reads
#optional: --rf for strandness

 kallisto quant\
 -i $WORKING_DIR/kallisto_genome_index\
 -b 50\
 -o $WORKING_DIR/H_1_1\
 --rf-stranded \
 $SAMPLE_FA/trim_1_1_L3_R1_001_ij43KLkHk1vK.fastq.gz $SAMPLE_FA/trim_1_1_L3_R2_001_qyjToP2TB6N7.fastq.gz

 kallisto quant\
 -i $WORKING_DIR/kallisto_genome_index\
 -b 50\
 -o $WORKING_DIR/H_1_2\
 --rf-stranded \
 $SAMPLE_FA/trim_1_2_L3_R1_001_DnNWKUYhfc9S.fastq.gz $SAMPLE_FA/trim_1_2_L3_R2_001_SNLaVsTQ6pwl.fastq.gz

 kallisto quant\
 -i $WORKING_DIR/kallisto_genome_index\
 -b 50\
 -o $WORKING_DIR/H_1_5\
 --rf-stranded \
 $SAMPLE_FA/trim_1_5_L3_R1_001_iXvvRzwmFxF3.fastq.gz $SAMPLE_FA/trim_1_5_L3_R2_001_iXCMrktKyEh0.fastq.gz

 kallisto quant\
 -i $WORKING_DIR/kallisto_genome_index\
 -b 50\
 -o $WORKING_DIR/P1\
 --rf-stranded \
 $SAMPLE_FA/trim_P1_L3_R1_001_9L0tZ86sF4p8.fastq.gz $SAMPLE_FA/trim_P1_L3_R2_001_yd9NfV9WdvvL.fastq.gz

 kallisto quant\
 -i $WORKING_DIR/kallisto_genome_index\
 -b 50\
 -o $WORKING_DIR/P2\
 --rf-stranded \
 $SAMPLE_FA/trim_P2_L3_R1_001_R82RphLQ2938.fastq.gz $SAMPLE_FA/trim_P2_L3_R2_001_06FRMIIGwpH6.fastq.gz

 kallisto quant\
 -i $WORKING_DIR/kallisto_genome_index\
 -b 50\
 -o $WORKING_DIR/P3\
 --rf-stranded \
 $SAMPLE_FA/trim_P3_L3_R1_001_fjv6hlbFgCST.fastq.gz $SAMPLE_FA/trim_P3_L3_R2_001_xo7RBLLYYqeu.fastq.gz
