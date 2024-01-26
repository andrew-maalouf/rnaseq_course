#!/usr/bin/Rscript

load("/data/courses/rnaseq_course/lncRNAs/Project1/references/Human_logitModel.RData")
test <- read.table(file="/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/coding_potential.dat",sep="\t",col.names=c("ID","mRNA","ORF","Fickett","Hexamer"))
test$prob <- predict(mylogit,newdata=test,type="response")
attach(test)
output <- cbind("mRNA_size"=mRNA,"ORF_size"=ORF,"Fickett_score"=Fickett,"Hexamer_score"=Hexamer,"coding_prob"=test$prob)
write.table(output,file="/data/courses/rnaseq_course/lncRNAs/Project1/users/amaalouf/step6/coding_potential.tsv",quote=F,sep="\t",row.names=ID)
