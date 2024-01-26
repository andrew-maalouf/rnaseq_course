#!/bin/bash

#SBATCH --mail-user=andrew.maalouf@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="QC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --partition=pall
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/error_%j.e

module load UHTS/Quality_control/fastqc/0.11.9; 
module load UHTS/Analysis/trimmomatic/0.36;

#create and go to the TP directory
mkdir /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/fastqc/PE/
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/fastqc/PE

ln -s /data/courses/rnaseq_course/lncRNAs/fastq/$1 .
ln -s /data/courses/rnaseq_course/lncRNAs/fastq/$2 .


fastqc -t 16 $1
fastqc -t 16 $2

trimmomatic PE -phred33 -threads 4 $1 $2 trim_$1 unpaired_$1 trim_$2 unpaired_$2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100

fastqc -t 16 trim_$1
fastqc -t 16 trim_$2


