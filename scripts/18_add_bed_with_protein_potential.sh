#!/usr/bin/env bash

#SBATCH --job-name="table"
#SBATCH --cpus-per-task=12
#SBATCH --time=33:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=pall

#get the coding potential column from coding_potential.tsv and add it to the full bed file
paste full_bed_file.bed <(cut -f 6 coding_potential.tsv) > bed_with_coding_potential.bed

#get coding potential only for novel genes. in other words extract lines from full bed file where transcript id is MSTRG
grep "MSTRG" bed_with_coding_potential.bed > novel_with_coding_potential.bed


#keep novels that are intergenic
awk 'NR==FNR{a[$2]=$8; next} $2 in a {print $0, a[$2]}' novel_with_coding_potential.bed intergenic_novel.bed > intergenic_novel_with_coding_potential.bed




#get intergenic novel genes that are non coding (coding potential < 0.364)
awk '$8 < 0.364' intergenic_novel_with_coding_potential.bed > intergenic_novel_non_coding.bed


#get intergenic novel genes that are coding (coding potential > 0.364)
awk '$8 > 0.364' intergenic_novel_with_coding_potential.bed > intergenic_novel_coding.bed


#get novel genes that are coding (codingpotential > 0.364)
awk '$8 > 0.364' novel_with_coding_potential.bed > novel_coding.bed




