#!/usr/bin/env Rscript

# Convert raw counts to normalization data.
# Requires NOISeq, edgeR
# Replace 0 with minius digit??(k = 1)
# Length correct is enabled.
# Usage: ./normalize count outprefix [read length = 150]
# Arguments:
#   infile: Raw count file from htseq-count or verse.
#   outprefix: Prefix of outfile.

output <- function(data,outfile,mysep="\t"){
  write.table(data,outfile,sep=mysep,quote=F)
}

Args <- commandArgs(TRUE)

# Read args
outprefix <- Args[2]
# Read length
if (length(Args) == 3){
  rl = as.integer(Args[3])
}else {
  rl = 150
}

# Read count file
mycounts0 <- read.table(Args[1],sep="\t",header=T)
row.names(mycounts0) <- mycounts0[,1]
mycounts0 <- mycounts0[,-1]

# Split to length and count file.
mycounts <- mycounts0[,1:(length(mycounts0[1,])-1)]
mylength <- mycounts0[,length(mycounts0[1,])]
names(mylength) <- row.names(mycounts)

# Using NOISeq do normalization.
library(NOISeq)
myRPKM <- rpkm(mycounts, long = mylength, k = NULL, lc = 1)
myUQUA <- uqua(mycounts, long = mylength, k = NULL, lc = 1)
myTMM <- tmm(mycounts, long = mylength, k = NULL, lc = 1)


library(edgeR)
# DESeq normalization.
rle.f = calcNormFactors(mycounts,method="RLE")
total <- colSums(mycounts)
rle.f = rle.f * (total/mean(total))
myRLE = t(t(mycounts)/rle.f)

# TPM
# Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. Günter P. Wagner, Koryu Kin and Vincent J. Lynch. Theory in Biosciences. 2012, 131(4): 281-285.
rl
total = colSums(mycounts * rl / mylength)
myTPM = mycounts * rl * 10**6 / (array(mylength,dim=c(length(mylength),1)) %*% array(total,dim=c(1,length(total))))

# Output normalize
output(myRPKM,paste(outprefix,"_rpkm.txt",sep=""))
output(myUQUA,paste(outprefix,"_uqua.txt",sep=""))
output(myTMM,paste(outprefix,"_tmm.txt",sep=""))
output(myRLE,paste(outprefix,"_rle.txt",sep=""))
output(myTPM,paste(outprefix,"_tpm.txt",sep=""))
