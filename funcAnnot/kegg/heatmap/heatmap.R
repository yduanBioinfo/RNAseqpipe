#!/usr/bin/env Rscript
inFile <- file("stdin",'r')
outFile <- commandArgs(TRUE)
mydata <- read.table(inFile,header=TRUE)
mydata <- as.matrix(mydata)

delzero <- function(data){
	mybool=(data[,1]==0);
	for(i in 2:length(data[1,])){
		mybool=mybool&(data[,i]==0)
	};
	return(data[!mybool,])
}

mydata <- delzero(mydata)

library(gplots)
pdf(paste(outFile,"pdf",sep="."))
heatmap.2(mydata,col=greenred(75),scale="row",density.info='none',trace='none',key.xlab='donw-up regulate')
dev.off()
