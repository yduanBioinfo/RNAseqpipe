#!/usr/bin/env Rscript

Usage="cat genelist | goplot.R - outdir ifonlylist
    genelist	can be pure list of gene ID or NOISeq out file
    genelist:
    GOBPID  GOCCID  GOMFID  Pvalue  OddsRatio       ExpCount        Count   Size    Term    GODomain        Direction       left    right   event_type
    GO:0000209                      0.009489876710779171    418.0   0.00950118764845606     1       2       protein polyubiquitination      BP      over    T4    T5      All
    GO:0070936                      0.009489876710779171    418.0   0.00950118764845606     1       2       protein K48-linked ubiquitination       BP      over    T4      T5      All
    GO:0006414                      0.0142178486596539      208.5   0.0142517814726841      1       3       translational elongation        BP      over    T4    T5      All
        GO:0005765              0.0169488797321661      174.5   0.0169971671388102      1       3       lysosomal membrane      CC      over    T4      T5      All

"

library(argparser)
library(ggplot2)

p <- arg_parser("Plot bubble for GO enrichment")
# Add a positional argument
p <- add_argument(p, "name",  help="Name of input file")
p <- add_argument(p, "--outfile", help="Name of output file", default="goplot.pdf")
p <- add_argument(p, "--width", help="width of graph", default=7)
p <- add_argument(p, "--height", help="height of graph", default=10)
p <- add_argument(p, "--units", help="unit of graph", default="in")
argv <- parse_args(p)

data <- read.table(ifelse(argv$name=="-","stdin",argv$name),header=TRUE,sep="\t")
# keep data$Term order
data$Term <- factor(data$Term, levels = rev(data$Term))
outfile <- argv$outfile

p<- ggplot(data, aes(GODomain, Term)) +
    geom_point(aes(size=Count, color=-1*log10(Pvalue))) +
    scale_color_gradient(low="orange",high="red",na.value="NA")+
    labs(color=expression(-log[10](Pvalue)),size="Count")+
    theme_bw()+
    theme(axis.text.y = element_text(size = rel(1.5)))
ggsave(outfile,width=argv$width,height=argv$height,units=argv$units)

