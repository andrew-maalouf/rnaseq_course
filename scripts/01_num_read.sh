#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4g
#SBATCH --time=00:15:00
#SBATCH --job-name=num_reads
#SBATCH --mail-user=andrew.maalouf@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/num_reads/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/num_reads/error_%j.e
#SBATCH --partition=pall


READS_DIR=/data/courses/rnaseq_course/lncRNAs/fastq/

count=$(echo $(zcat $READS_DIR/$1|wc -l)/4|bc)

echo "The number of reads in $1 is: $count"
