#!/usr/bin/env bash

#SBATCH --job-name="stat"
#SBATCH --partition=pall
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/my_stat.o

READ_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/merge_assembly/stringtie_merged.gtf

exons_total=$(awk '$3 == "exon" {print}' $READ_DIR | wc -l)
#or faster: awk '$3 == "exon" {count++} END {print "Total number of exons:", count}' merge_assembly/stringtie_merged.gtf

transcripts_total=$(awk '$3 == "transcript" {print}' $READ_DIR | wc -l)
#or faster: awk '$3 == "transcript" {count++} END {print "Total number of transcripts:", count}' merge_assembly/stringtie_merged.gtf

genes_total=$(awk '$3 == "transcript" {print $10}' $READ_DIR | sort | uniq -c | wc -l)


single_exon_transcript=$(awk '$3 == "exon" {print $12}' $READ_DIR | sort | uniq -c | awk '$1 == 1' | wc -l)
#or : awk '$3 == "exon" {transcript_id = $12; exon_count[transcript_id]++} END {for (id in exon_count) {if (exon_count[id] == 1) {print id}}}' merge_assembly/stringtie_merged.gtf | wc -l

single_transcript_gene=$(awk '$3 == "transcript" {print $10}' $READ_DIR | sort | uniq -c | awk '$1 == 1' | wc -l)

single_exon_gene=$(awk '$3 == "exon" {gene_id = $10; exon_count[gene_id]++} END {for (id in exon_count) {if (exon_count[id] == 1) {print id}}}' $READ_DIR | wc -l)

unique_novel_exons=$(awk ' $3 == "exon" {if($12 !~ /ENST/) print $10}' $READ_DIR | sort | uniq -c | awk '$1 == 1 {print $2}' | wc -l)
novel_exons=$(awk '$12 !~ /ENST/ && $3 == "exon" {print $10}' $READ_DIR | wc -l)

unique_novel_transcripts=$(awk '$12 !~ /ENST/ && $3 == "exon" {print $12}' $READ_DIR | sort | uniq -c | awk '$1 == 1 {print $2}' | wc -l)
novel_transcripts=$(awk '$3 == "transcript" {if($12 !~ /ENST/) print $12}' $READ_DIR | sort | uniq | wc -l)

novel_genes=$(awk '$12 !~ /ENST/ && $3 == "transcript" {print $10}' $READ_DIR | sort | uniq | wc -l)

echo "The total number of exons is $exons_total"
echo "The total number of transcripts is $transcripts_total"
echo "The total number of genes is $genes_total"
echo "The total number of single exon transcript is $single_exon_transcript"
echo "The total number of single transcript gene is $single_transcript_gene"
echo "The total number of single exon gene is $single_exon_gene"
echo "The total number of novel genes is $novel_genes"
echo "The total number of novel exons is $novel_exons, of which $unique_novel_exons are unique (belong to single-exon gene)"
echo "The total number of novel transcripts is $novel_transcripts, of which $unique_novel_transcripts are single-exons"

awk '$12 !~ /ENST/ && $3 == "exon" {print $12}' $READ_DIR | sort | uniq -c | awk '$1 > 1 {print $0}' > novel_transcripts_with_multiple_exons.txt








echo ""
echo ""
echo "Now considering annotated chromosomes only:"
echo ""
echo ""

exons_total=$(awk '$3 == "exon" && $1 ~ /chr/ {print}' $READ_DIR | wc -l)
#or faster: awk '$3 == "exon" {count++} END {print "Total number of exons:", count}' merge_assembly/stringtie_merged.gtf

transcripts_total=$(awk '$3 == "transcript" && $1 ~ /chr/ {print}' $READ_DIR | wc -l)
#or faster: awk '$3 == "transcript" {count++} END {print "Total number of transcripts:", count}' merge_assembly/stringtie_merged.gtf

genes_total=$(awk '$3 == "transcript" && $1 ~ /chr/ {print $10}' $READ_DIR | sort | uniq -c | wc -l)


single_exon_transcript=$(awk '$3 == "exon" && $1 ~ /chr/ {print $12}' $READ_DIR | sort | uniq -c | awk '$1 == 1' | wc -l)
#or : awk '$3 == "exon" {transcript_id = $12; exon_count[transcript_id]++} END {for (id in exon_count) {if (exon_count[id] == 1) {print id}}}' merge_assembly/stringtie_merged.gtf | wc -l

single_transcript_gene=$(awk '$3 == "transcript" && $1 ~ /chr/ {print $10}' $READ_DIR | sort | uniq -c | awk '$1 == 1' | wc -l)

single_exon_gene=$(awk '$3 == "exon" && $1 ~ /chr/ {gene_id = $10; exon_count[gene_id]++} END {for (id in exon_count) {if (exon_count[id] == 1) {print id}}}' $READ_DIR | wc -l)

unique_novel_exons=$(awk '$12 !~ /ENST/ && $3 == "exon" && $1 ~ /chr/ {print $10}' $READ_DIR | sort | uniq -c | awk '$1 == 1 {print $2}' | wc -l)
novel_exons=$(awk '$12 !~ /ENST/ && $3 == "exon" && $1 ~ /chr/ {print $10}' $READ_DIR | wc -l)

unique_novel_transcripts=$(awk '$12 !~ /ENST/ && $3 == "exon" && $1 ~ /chr/ {print $12}' $READ_DIR | sort | uniq -c | awk '$1 == 1 {print $2}' | wc -l)
novel_transcripts=$(awk '$12 !~ /ENST/ && $3 == "transcript" && $1 ~ /chr/ {print $12}' $READ_DIR | sort | uniq | wc -l)

novel_genes=$(awk '$12 !~ /ENST/ && $3 == "transcript" && $1 ~ /chr/ {print $10}' $READ_DIR | sort | uniq | wc -l)

echo "The total number of exons is $exons_total"
echo "The total number of transcripts is $transcripts_total"
echo "The total number of genes is $genes_total"
echo "The total number of single exon transcript is $single_exon_transcript"
echo "The total number of single transcript gene is $single_transcript_gene"
echo "The total number of single exon gene is $single_exon_gene"
echo "The total number of novel genes is $novel_genes"
echo "The total number of novel exons is $novel_exons, of which $unique_novel_exons are unique (belong to single-exon gene)"
echo "The total number of novel transcripts is $novel_transcripts, of which $unique_novel_transcripts are single-exons"
