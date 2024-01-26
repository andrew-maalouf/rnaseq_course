#!/usr/bin/env bash

#SBATCH --mail-user=andrew.maalouf@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="index_genome"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/index_output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/index_error_%j.e

THREADS=$SLURM_CPUS_PER_TASK

module load UHTS/Aligner/hisat/2.2.1;

REFERENCE_READ=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/genome/ 
WORK_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2

mkdir $WORK_DIR/hisat_reference_index
cd $WORK_DIR/hisat_reference_index

hisat2-build -p $THREADS $REFERENCE_READ/GRCh38.p14.genome.fa $WORK_DIR/hisat_reference_index
