#!/usr/bin/env Rscript

help = "Usage:./mykeggenrich.R infile\n Outfile will be in the same directory"
File <- commandArgs(TRUE)
tmp <- tryCatch(Table.a <- read.table(file=File,head=TRUE,sep="\t"), error=function(e) NULL)
if(is.null(tmp)){print("Warning: The input file is empty");q()}
fisherMul <- function(X1){
dim(X1)<-c(2,2)
X1_fisher <- fisher.test(X1)
X1_fisher[[1]]
}

#p.value<-sapply(data.frame(rbind(Table.a[,2],Table.a[,3],Table.a[,4],Table.a[,5])),FUN=fisherMul)
p.value<-sapply(data.frame(rbind(Table.a[,2],(Table.a[,3]-Table.a[,2]),Table.a[,4],(Table.a[,5]-Table.a[,4]))),FUN=fisherMul)
Fisher.value<-as.vector(p.value)
pvalue.adj <- as.vector(p.adjust(Fisher.value,method="fdr"))#FDR adjust
file.out<-data.frame(Table.a[,1:5],Fisher.value,pvalue.adj)
Path.Diff<-file.out[file.out[,6]<0.05,]
Path.adjDiff<-file.out[file.out[,7]<0.05,]
write.table(file.out,file=paste(File,"OUT",sep="."),sep="\t",quote=FALSE,row.names=FALSE)
write.table(Path.Diff,file=paste(File,"Diff",sep="."),sep="\t",quote=FALSE,row.names=FALSE)
write.table(Path.adjDiff,file=paste(File,"adjDiff",sep="."),sep="\t",quote=FALSE,row.names=FALSE)
