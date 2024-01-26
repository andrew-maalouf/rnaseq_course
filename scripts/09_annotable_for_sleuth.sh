#!/usr/bin/env bash

#SBATCH --job-name="table"
#SBATCH --partition=pall


MERGED_ASSEMBLY=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/merge_assembly/stringtie_merged.gtf
GENCODE_ANNOTATION=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/genome/gencode.v44.annotation.gtf
WORK_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3

# Extract gene_id, transcript_id, and gene_name from MERGED_ASSEMBLY, sort, and remove quotes
#if name is ot available, write "Unknown"
awk '$3=="transcript" {if($15=="ref_gene_id") print $12,$16,$14; else print $12,$10,$14;}' $MERGED_ASSEMBLY | sort | tr -d "'" | awk '{if ($3 == "") $3 = "unknown"; print}'> $WORK_DIR/table_sleuth.txt

#remove ; and ""
sed -i 's/[";]//g' table_sleuth.txt
