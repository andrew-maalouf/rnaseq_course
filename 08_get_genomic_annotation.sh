#!/usr/bin/env bash

#SBATCH --job-name="table"
#SBATCH --partition=pall


MERGED_ASSEMBLY=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/merge_assembly/stringtie_merged.gtf
GENCODE_ANNOTATION=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/genome/gencode.v44.annotation.gtf
WORK_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3

# Extract gene_id, transcript_id, and gene_name from MERGED_ASSEMBLY, sort, and remove quotes
awk '$3=="transcript" {if($15=="ref_gene_id") print $12,$16,$14; else print $12,$10,$14;}' $MERGED_ASSEMBLY | sort | tr -d "'" > $WORK_DIR/table_annotation.txt

# Extract biotype per transcript_id from GENCODE_ANNOTATION, sort, and remove quotes
awk ' $3=="transcript"  {print $12,$14,$18}' $GENCODE_ANNOTATION | sort | tr -d "'" > $WORK_DIR/transcript_biotype.txt

# Join the two tables, remove quotes, and add "novel" for missing columns
join -t ';' -a 1 -a 2 -1 1 -2 1 $WORK_DIR/table_annotation.txt $WORK_DIR/transcript_biotype.txt | sed 's/;;/;/g' | awk '{gsub(/'\''/,"",$3); gsub(/'\''/,"",$4); print $1, $2, $3, $4}' | awk '{if ($4 == "") $4 = "novel"; print}' > $WORK_DIR/final_table.txt

# Keep lines that contain "novel", "protein_coding", "lncRNA"
grep -E 'novel|protein_coding|lncRNA' final_table.txt > filtered_final_table.txt

# Delete residual files
rm transcript_biotype.txt
rm table_annotation.txt
rm final_table.txt
