#!/usr/bin/env bash

#SBATCH --mail-user=andrew.maalouf@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="merge_assembly"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/index_output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/index_error_%j.e

WORK_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/merge_assembly
READ_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/transcriptome_assembly
REF_ANNOT=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/genome/gencode.v44.annotation.gtf


#load modules
module load UHTS/Aligner/stringtie/1.3.3b;
module load UHTS/Assembler/cufflinks/2.2.1;

#create output directory
mkdir merge_assembly

#enter directory where gtf files are found instead of creating links in output directory
cd $READ_DIR

#the following step has been included in the previous script
#assembly_GTF_list.txt is a text file “manifest” with a list (one per line) of GTF files that I would like to merge together into a single GTF file.
#for file in `ls -1 *.gtf`; do echo $file >> assembly_GTF_list.txt; done

#run stringtie
#-–rf tells StringTie that our data is stranded and to use the correct strand specific mode (i.e. assume a stranded library fr-firststrand).
#-p 6 tells stringtie to use six CPUs
#-o tells stringtie to write output to a particular file or directory
#-G tells stringtie where to find reference gene annotations. It will use these annotations to gracefully merge novel isoforms (for de novo runs) and known isoforms and maximize overall assembly quality.

stringtie --rf --merge -p 6 -o $WORK_DIR/stringtie_merged.gtf -G $REF_ANNOT $READ_DIR/assembly_GTF_list.txt


#Compare reference guided transcripts to the known annotations. This allows us to assess the quality of transcript predictions made from assembling the RNA-seq data.
cd $WORK_DIR
cuffcompare -r $REF_ANNOT -o gffcompare stringtie_merged.gtf
