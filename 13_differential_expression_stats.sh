#!/usr/bin/env bash
output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step5/differential_expression_stats.txt
TRANSCRIPT=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step5/transcript_level_results.tsv

# note: when looking for "novel", the column before is empty,
# so the $ of column will be biotype-1 and conditions -1

echo "Disclaimer: upregulated means upregulated in holoclones./// differentially expressed initial conditions are qval < 0.05 and |log2fc| > 1" > $output
echo "" >> $output

echo "Transcript level stats" >> $output
all=$(cat $TRANSCRIPT | wc -l)
tested=$(grep -v "NA" $TRANSCRIPT | wc -l)
signif=$(awk '{if ($3 == "novel") print ($5 < 0.05 && $6 >= 1 || $5 < 0.05 && $6 <= -1) ; else print (($6 < 0.05 && $7 >= 1 || $6 < 0.05 && $7 <= -1) && ($4 == "lncRNA" || $4 == "protein_coding"))}' $TRANSCRIPT | awk '$1' | wc -l)


diff_lncRNAs=$(awk '$4 == "lncRNA" {print}' $TRANSCRIPT | awk '$6 < 0.05 && $7 >= 1 || $6 < 0.05 && $7 <= -1 {print}' | wc -l)
up_r=$(awk '$4 == "lncRNA" {print}' $TRANSCRIPT | awk '$6 < 0.05 && $7 <= -1 {print}' | wc -l)
down_r=$(awk '$4 == "lncRNA" {print}' $TRANSCRIPT | awk '$6 < 0.05 && $7 >= 1 {print}' | wc -l)

diff_protein=$(awk '$4 == "protein_coding" {print}' $TRANSCRIPT | awk '$6 < 0.05 && $7 >= 1 || $6 < 0.05 && $7 <= -1 {print}' | wc -l)
up_p=$(awk '$4 == "protein_coding" {print}' $TRANSCRIPT | awk '$6 < 0.05 && $7 <= -1 {print}' | wc -l)
down_p=$(awk '$4 == "protein_coding" {print}' $TRANSCRIPT | awk '$6 < 0.05 && $7 >= 1 {print}' | wc -l)

diff_novel=$(awk '$3 == "novel" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 >= 1 || $5 < 0.05 && $6 <= -1 {print}' | wc -l)
up_n=$(awk '$3 == "novel" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 <= -1 {print}' | wc -l)
down_n=$(awk '$3 == "novel" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 >= 1 {print}' | wc -l)

echo "Out of $all transcripts (-1), $tested have values for differential expression" >> $output
echo "$signif are differentially expressed of which:" >> $output
echo "Out of all differentially expressed transcripts:" >> $output
echo "$diff_lncRNAs are lncRNAs, upregulated: $up_r & downregulated: $down_r" >> $output
echo "$diff_protein are protein coding, upregulated: $up_p & downregulated: $down_p" >> $output
echo "$diff_novel are novel, upregulated: $up_n & downregulated: $down_n" >> $output

# I am lazy to change "transcript" to "gene"
TRANSCRIPT=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step5/gene_level_results.txt

echo "" >> $output
echo "Gene Level" >> $output
all=$(cat $TRANSCRIPT | wc -l)
tested=$(grep -v "NA" $TRANSCRIPT | wc -l)
signif=$(awk '{if ($2 == "novel") print ($4 < 0.05 && $5 >= 1 || $4 < 0.05 && $5 <= -1) ; else print ($5 < 0.05 && $6 >= 1 || $5 < 0.05 && $6 <= -1)}' $TRANSCRIPT | awk '$1' | wc -l)


diff_lncRNAs=$(awk '$3 == "lncRNA" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 >= 1 || $5 < 0.05 && $6 <= -1 {print}' | wc -l)
up_r=$(awk '$3 == "lncRNA" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 <= -1 {print}' | wc -l)
down_r=$(awk '$3 == "lncRNA" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 >= 1 {print}' | wc -l)

diff_protein=$(awk '$3 == "protein_coding" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 >= 1 || $5 < 0.05 && $6 <= -1 {print}' | wc -l)
up_p=$(awk '$3 == "protein_coding" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 <= -1 {print}' | wc -l)
down_p=$(awk '$3 == "protein_coding" {print}' $TRANSCRIPT | awk '$5 < 0.05 && $6 >= 1 {print}' | wc -l)

diff_novel=$(awk '$2 == "novel" {print}' $TRANSCRIPT | awk '$4 < 0.05 && $5 >= 1 || $4 < 0.05 && $5 <= -1 {print}' | wc -l)
up_n=$(awk '$2 == "novel" {print}' $TRANSCRIPT | awk '$4 < 0.05 && $5 <= -1 {print}' | wc -l)
down_n=$(awk '$2 == "novel" {print}' $TRANSCRIPT | awk '$4 < 0.05 && $5 >= 1 {print}' | wc -l)


echo "Out of $all genes (-1), $tested have values for differential expression" >> $output
echo "$signif are differentially expressed of which:" >> $output
echo "Out of all differentially expressed genes:" >> $output
echo "$diff_lncRNAs are lncRNAs, upregulated: $up_r & downregulated: $down_r" >> $output
echo "$diff_protein are protein coding, upregulated: $up_p & downregulated: $down_p" >> $output
echo "$diff_novel are novel, upregulated: $up_n & downregulated: $down_n" >> $output
