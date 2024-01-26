📝 In this study, I aim to identify protein-coding genes and lncRNAs that phenotypically distinguishes the holoclones in the A549 cell line. Therefore, following RNA-seq experiment producing reads, the latter will be used to annotate transcripts from known and novel genes. Next, this annotation will be used to quantify and identify differentially expressed genes and transcripts between the holoclones and the parental line. I aspire to discover candidates that could be targeted to alter the phenotypic trajectories of NSCLC cells.

🔗 Link to fastqc files:
  /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step1/fastqc/PE

🔗 Link to bam files:
  /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step2/read_align

🔗 Link to merged transcriptome assembly:
  /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step3/merge_assembly/stringtie_merged.gtf

🔗 Link to kallisto results:
  /data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step4/kallisto_expression

📂 Files content:
  1- abundance_stats.txt: quantification results statistics
  2- differential_expression_stats.txt: differential expression results statistics
  3- end_table.tsv: table including all 20066 novel transcripts on reference chromosomes with their respective transcript level differential expression results + 5' 3' completeness + intergenic + coding potential
  4- gene_level_results.txt: differential expression full sleuth results on gene level
  5- lnc_diff_signif_intergenic.tsv: list of all intergenic novel transcripts with qval < 0.05 filtered out from end_table.tsv
  6- my_stat.o: statistics of merged transcriptome assembly
  7- possible_candidates.tsv: list of all intergenic novel transcripts with qval < 0.05 and |log2fc| > 1 filtered out from end_table
  8- transcript_level_results.tsv: differential expression full sleuth results on transcript level
  9- written_answers.txt: explicit answers to step6 questions


