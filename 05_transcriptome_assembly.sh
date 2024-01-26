#!/usr/bin/env bash

#SBATCH --mail-user=andrew.maalouf@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="transcriptome_assembly"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/index_output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/index_error_%j.e

WORK_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/transcriptome_assembly
READ_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/read_align/
REF_ANNOT=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/genome/gencode.v44.annotation.gtf
BASE_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/

#load modules
module load UHTS/Aligner/stringtie/1.3.3b;

#create and enter directory
mkdir $BASE_DIR/transcriptome_assembly
cd $WORK_DIR

#create bam files links into my directory
ln -s $READ_DIR/*_sorted.bam .

#create gtf file for each bam file through a loop
#--rf : tells StringTie that our data is stranded and to use the correct strand specific mode
#-p 4 tells Stringtie to use eight CPUs
#-G reference annotation to use for guiding the assembly process (GTF/GFF3)
#-l name prefix for output transcripts (default: STRG)
#-o output path/file name for the assembled transcripts GTF (default: stdout)


for file in `ls -1 *.bam`;
do bam="${file%_sorted.bam}";
echo "working on $bam";
stringtie --rf -p 6 -G $REF_ANNOT -l $bam -o $WORK_DIR/${bam}_transcript.gtf $file;
done

cd $WORK_DIR

#assembly_GTF_list.txt is a text file “manifest” with a list (one per line) of GTF files that I would like to merge together into a single GTF file.
for file in `ls -1 *.gtf`;
do echo $file >> assembly_GTF_list.txt;
done

