#!/usr/bin/env Rscript-3.6.1

Usage = "./format_godb.R xxxx.godb [species_name]"
Help = "Save godb as GeneSetCollection(gsc) object."
library(GOstats)
library("GSEABase")

Args <- commandArgs(TRUE)
dbstatsFile = Args[1]
species <- ifelse(is.na(Args[2]), "grass_carp", Args[2])

outdbfile=paste(dbstatsFile,'RData',sep='.')
universeData = read.table(dbstatsFile,header=F,col.names=c("go_id","evidence","gene_id"))
goFrame = GOFrame(universeData,organism="grass_carp")
goAllFrame = GOAllFrame(goFrame)
gsc = GeneSetCollection(goAllFrame,setType=GOCollection())
save(gsc, file=outdbfile)
