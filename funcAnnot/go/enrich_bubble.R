#!/usr/bin/env Rscript-4.0.2

Args <- commandArgs(TRUE)
#two arguments: 
#go_annot_db_file(gostats format) DE_out_file(NOISeq out format)
#test gene in godb
#version: 0.0.1
Usage="cat genelist | goplot.R - outdir ifonlylist
    genelist	can be pure list of gene ID or NOISeq out file
    genelist:
    GOBPID  GOCCID  GOMFID  Pvalue  OddsRatio       ExpCount        Count   Size    Term    GODomain        Direction       left    right   event_type
    GO:0000209                      0.009489876710779171    418.0   0.00950118764845606     1       2       protein polyubiquitination      BP      over    T4    T5      All
    GO:0070936                      0.009489876710779171    418.0   0.00950118764845606     1       2       protein K48-linked ubiquitination       BP      over    T4      T5      All
    GO:0006414                      0.0142178486596539      208.5   0.0142517814726841      1       3       translational elongation        BP      over    T4    T5      All
        GO:0005765              0.0169488797321661      174.5   0.0169971671388102      1       3       lysosomal membrane      CC      over    T4      T5      All

"

data <- read.table(ifelse(Args[1]=="-","stdin",Args[1]),header=TRUE,sep="\t")
# keep data$Term order
data$Term <- factor(data$Term, levels = rev(data$Term))
outfile <- ifelse(is.na(Args[2]),"goplot.pdf",Args[2])

library(ggplot2)
p<- ggplot(data, aes(GODomain, Term)) +
    geom_point(aes(size=Count, color=-1*log10(Pvalue))) +
    scale_color_gradient(low="orange",high="red",na.value="NA")+
    labs(color=expression(-log[10](Pvalue)),size="Count")+
    theme_bw()+
    theme(axis.text.y = element_text(size = rel(1.5)))
ggsave(outfile)

