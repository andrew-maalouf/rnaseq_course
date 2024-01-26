#!/usr/bin/env bash

#SBATCH --mail-user=andrew.maalouf@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="align_reads"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/index_output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/index_error_%j.e

XX=$1
THREADS=$SLURM_CPUS_PER_TASK

#index genome reference before running this script

# load modules for alignment
module add UHTS/Aligner/hisat/2.2.1;
module load UHTS/Analysis/samtools/1.10;

READ_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/fastqc/PE
INDEX_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/hisat_reference_index
WORK_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2

mkdir read_align
cd read_align

ln -s $READ_DIR/trim_${XX}_*_R1_*.fastq.gz .
ln -s $READ_DIR/trim_${XX}_*_R2_*.fastq.gz .


hisat2 -p $THREADS \
 --rna-strandness R \
 -x $INDEX_DIR/hisat_reference_index \
  -1 trim_${XX}_*_R1_*.fastq.gz \
   -2 trim_${XX}_*_R2_*.fastq.gz \
    -S ${XX}.sam

#transform sam to bam
samtools view -b -@ $THREADS ${XX}.sam > ${XX}_unsorted.bam

#sort bam
samtools sort -@ $THREADS ${XX}_unsorted.bam > ${XX}_sorted.bam

#index bam > bai
samtools index ${XX}_sorted.bam

#cleaning directory
cd $WORK_DIR/read_align
rm *.sam
rm *_unsorted.bam