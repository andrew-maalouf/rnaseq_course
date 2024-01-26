#!/usr/bin/env bash
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/quality_output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/quality_error_%j.e

WORKING_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_expression
BASE_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/

cd $WORKING_DIR
#sum tpm values and divide them by million
#expect to get value of 1 if quality is good
#create loop to go through all of the 6 abundance.tsv files at once

for each_replicate_dir in "$WORKING_DIR"/*/; do
    if [ -f "$each_replicate_dir/abundance.tsv" ]; then
        sum_tpm=$(awk '{ sum += $5 } END { print sum / 1000000 }' "$each_replicate_dir/abundance.tsv" );
        #num_transcript=$(awk 'NR>1 {print $1}' "$each_replicate_dir/abundance.tsv" )
        echo "The tpm sum for $each_replicate_dir is $sum_tpm";
        cd $WORKING_DIR;
        fi;
    done