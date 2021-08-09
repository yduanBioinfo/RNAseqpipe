#!/usr/bin/env Rscript-3.6.1

#####!/usr/bin/env Rscript-3.2.5

GOstat_and_output <- function(genelist,dbstatsFile,outdir)
{
    library(GOstats)
    library("GSEABase")
    #get gene id
    universeData = read.table(dbstatsFile,header=F,col.names=c("go_id","evidence","gene_id"))
    universID=as.character(unique(sort(universeData$gene_id)))
    geneID=unique(sort(intersect(universID,genelist$ID)))
    
    formated.godb = paste(dbstatsFile,'RData',sep='.')
    if (file.exists(formated.godb))
    {
        attach(formated.godb)
        }else{
        goFrame = GOFrame(universeData,organism="grass_carp")
        goAllFrame = GOAllFrame(goFrame)
        gsc = GeneSetCollection(goAllFrame,setType=GOCollection())
    }
    GOstat_one_file(geneID, universID, gsc, outdir)
}

GOstat_one_file<-function(geneID, universID, gsc, outdir)
{
    # gsc: GeneSetCollection

    params = GSEAGOHyperGParams(name="my_params",geneSetCollection=gsc,geneIds=geneID,universeGeneIds=universID,ontology="MF",pvalueCutoff=0.05,conditional=F,testDirection="over")
    #test and write out for cutoff = 0.05
    MF_over = hyperGTest(params)
    testDirection(params)<-"under"
    MF_under = hyperGTest(params)
    testDirection(params)<-"over"
    ontology(params)<-"BP"
    BP_over = hyperGTest(params)
    testDirection(params)<-"under"
    BP_under = hyperGTest(params)
    testDirection(params)<-"over"
    ontology(params)<-"CC"
    CC_over = hyperGTest(params)
    testDirection(params)<-"under"
    CC_under = hyperGTest(params)
    #write pvalue cutoff = 0.05
    write.table(summary(MF_over),file=paste(outdir,"MF_over.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(summary(MF_under),file=paste(outdir,"MF_under.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(summary(BP_over),file=paste(outdir,"BP_over.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(summary(BP_under),file=paste(outdir,"BP_under.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(summary(CC_over),file=paste(outdir,"CC_over.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(summary(CC_under),file=paste(outdir,"CC_under.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}

Args <- commandArgs(TRUE)
#two arguments: 
#go_annot_db_file(gostats format) DE_out_file(NOISeq out format)
#test gene in godb
#version: 0.0.1
Usage="cat genelist | gostat_goplot_RamiGO_out.R - dbfile outdir ifonlylist
    genelist	can be pure list of gene ID or NOISeq out file
    dbfile	
        GO:0005783	ISA	CIWT.295.1
        GO:0055074	ISA	CIWT.295.1
        GO:0005783	ISA	CIWT.295.2
        GO:0055074	ISA	CIWT.295.2
        GO:0005783	ISA	CIWT.295.3
        if dbfile.RData exist, use RData instead.
    outdir	outprefix (defalut:./goenrich)
    ifonlylist	type of genelist. T,True,[empty] means pure list of gene ID, others means NOISeq output.

    When the provided genelist has not enriched to any GO term, the fowllowing error will be reported:
        Error in order(pvals) : argument 1 is not a vector
        Calls: GOstat_and_output ... .hyperGTestInternal -> new -> initialize -> initialize -> order
        Execution halted
    This bug should be fixed in the future.
"

genelist <- read.table(ifelse(Args[1]=="-","stdin",Args[1]))
#read file
dbstatsFile <- Args[2]#go annot db file
#read genelist
outdir <- ifelse(is.na(Args[3]),"goenrich",Args[3])
outdir <- ifelse(("/"==substr(outdir,nchar(outdir),nchar(outdir))),outdir,paste(outdir,"/",sep=""))
dir.create(outdir)
# When apply NOISeq result, set ifonlylist to ???True(any word except for T or True).???
# When apply genelist, set ifonlylist to True or leave ifonlylist empty.
ifonlylist <- ifelse(is.na(Args[4]) | (Args[4] == "T") | (Args[4] == "True"),T,F)
if(ifonlylist){
    genelist$ID <- genelist[,1]
}else{
    genelist$ID <- rownames(genelist)
}

GOstat_and_output(genelist,dbstatsFile,outdir)
