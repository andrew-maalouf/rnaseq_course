#!/usr/bin/env bash
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/abundance_stat_output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/abundance_stat_error_%j.e

DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_expression
output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/abundance_stats.txt

cd $DIR

for sample in $(find . -maxdepth 1 -type d -printf "%P\n");
do

    echo "Stats for $sample" >> $output

    #total estimated counts
    counts=$(awk -F'\t' 'BEGIN{n=0} NR>1 {n=n+$4} END{print n}' $DIR/$sample/abundance.tsv)
    echo "total of estimated counts: $counts" >> $output

    #number of detected transcripts
    total_transcripts=$(awk -F'\t' 'NR>1 {if($4 > 0) ++n} END{print n}' $DIR/$sample/abundance.tsv)
    echo "number of transcripts detected: $total_transcripts" >> $output

    #number of detected genes
    total_genes=$(awk -F'\t' '$1 ~ /\.1$/ {if($4 > 0) ++n } END{print n}' $DIR/$sample/abundance.tsv)
    echo "number of genes detected: $total_genes" >> $output

    #number of detected novel transcripts 
    novel_transcripts=$(awk -F'\t' 'NR>1 && $1 !~ /'ENS'/ {if($4 > 0) ++n} END{print n}' $DIR/$sample/abundance.tsv) 
    echo "number of novel transcripts detected: $novel_transcripts" >> $output

    #number of detected novel genes 
    novel_genes=$(awk -F'\t' '$1 ~ /\.1$/ && $1 !~ /'ENS'/ {if($4 > 0) ++n } END{ print n}' $DIR/$sample/abundance.tsv)
    echo "number of novel genes detected: $novel_genes" >> $output;
done

