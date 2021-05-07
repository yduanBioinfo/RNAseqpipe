#!/usr/bin/env Rscript-3.2.5

GOstat_and_output <- function()
{
    library(GOstats)
    library("GSEABase")
    #get gene id
    universID=as.character(unique(sort(universeData$gene_id)))
    geneID=unique(sort(intersect(universID,genelist$ID)))
    
    #creat param obj
    goFrame = GOFrame(universeData,organism="grass_carp")
    goAllFrame = GOAllFrame(goFrame)
    gsc = GeneSetCollection(goAllFrame,setType=GOCollection())
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
    params
}
prepare_for_plot<-function(params)
{
    #set no pvalue cutoff to do fdr adjust
    pvalueCutoff(params) <- 1
    ontology(params)<-"BP"
    testDirection(params)<-"over"
    BP_over1 = hyperGTest(params)
    ontology(params)<-"MF"
    MF_over1 = hyperGTest(params)
    ontology(params)<-"CC"
    CC_over1 = hyperGTest(params)
    ontology(params)<-"BP"
    testDirection(params)<-"under"
    BP_under1 = hyperGTest(params)
    ontology(params)<-"MF"
    MF_under1 = hyperGTest(params)
    ontology(params)<-"CC"
    CC_under1 = hyperGTest(params)
    #get data frame
    gostatBPover <- summary(BP_over1)
    gostatMFover <- summary(MF_over1)
    gostatCCover <- summary(CC_over1)
    gostatBPunder <- summary(BP_under1)
    gostatMFunder <- summary(MF_under1)
    gostatCCunder <- summary(CC_under1)
    #fdr adjust and join into one data frame
    colnames(gostatBPover) <- c("ID","Pvalue","OddsRatio","ExpCount","Count","Size","term")
    colnames(gostatMFover) <- c("ID","Pvalue","OddsRatio","ExpCount","Count","Size","term")
    colnames(gostatCCover) <- c("ID","Pvalue","OddsRatio","ExpCount","Count","Size","term")
    colnames(gostatBPunder) <- c("ID","Pvalue","OddsRatio","ExpCount","Count","Size","term")
    colnames(gostatMFunder) <- c("ID","Pvalue","OddsRatio","ExpCount","Count","Size","term")
    colnames(gostatCCunder) <- c("ID","Pvalue","OddsRatio","ExpCount","Count","Size","term")
    gostatBPover$adj_pval <- p.adjust(gostatBPover[,2],method="fdr")
    gostatMFover$adj_pval <- p.adjust(gostatMFover[,2],method="fdr")
    gostatCCover$adj_pval <- p.adjust(gostatCCover[,2],method="fdr")
    gostatBPunder$adj_pval <- p.adjust(gostatBPunder[,2],method="fdr")
    gostatMFunder$adj_pval <- p.adjust(gostatMFunder[,2],method="fdr")
    gostatCCunder$adj_pval <- p.adjust(gostatCCunder[,2],method="fdr")
    gostatBPover$category <- rep("BP",length(gostatBPover[,2]))
    gostatBPunder$category <- rep("BP",length(gostatBPunder[,2]))
    gostatMFover$category <- rep("MF",length(gostatMFover[,2]))
    gostatMFunder$category <- rep("MF",length(gostatMFunder[,2]))
    gostatCCover$category <- rep("CC",length(gostatCCover[,2]))
    gostatCCunder$category <- rep("CC",length(gostatCCunder[,2]))
    mygoterms = rbind(gostatCCover,gostatCCunder,gostatMFover,gostatMFunder,gostatBPover,gostatBPunder)
    #function get_geneids2gostat
    #find the genes count to the same GO annotion
    get_offspring <- function(ID,category){if(category=="CC"){return(unlist(c(ID,get(ID,GOCCOFFSPRING))))};if(category=="BP"){return(unlist(c(ID,get(ID,GOBPOFFSPRING))))};if(category=="MF"){return(unlist(c(ID,get(ID,GOMFOFFSPRING))))}}
    get_geneid  <- function(goid,godb=goData){return(as.vector(godb[godb$go_id==goid,3]))}
    get_geneids <- function(go){ids=unlist(lapply(go,get_geneid));myids=ids[!is.na(ids)];return(myids)}
    get_geneids2gostat <- function(go,cate){gos=get_offspring(go,cate);genes=unlist(lapply(gos,get_geneids));return(paste(genes,collapse=","))}
    gos=get_offspring(gostatBPover$ID,"BP")
    genes=unlist(lapply(gos,get_geneids))
    #add genes 
    gostatBPover$genes = lapply(gostatBPover$ID,get_geneids2gostat,"BP")
    gostatBPunder$genes = lapply(gostatBPunder$ID,get_geneids2gostat,"BP")
    gostatMFover$genes = lapply(gostatMFover$ID,get_geneids2gostat,"MF")
    gostatMFunder$genes = lapply(gostatMFunder$ID,get_geneids2gostat,"MF")
    gostatCCover$genes = lapply(gostatCCover$ID,get_geneids2gostat,"CC")
    gostatCCunder$genes = lapply(gostatCCunder$ID,get_geneids2gostat,"CC")
    gostatall <- rbind(gostatBPover,gostatBPunder,gostatCCover,gostatCCunder,gostatMFover,gostatMFunder)
}
plots <-function()
{
    # circle_dat function requires NOISeq-like genelist, hance,
    # only on that circumstances, can we plot.

    library(GOplot)
    #sort zero genes to the last line
    mysub = subset(gostatall,adj_pval<0.05)
    if(length(mysub[,1])<1){
    	stop("no sig genes")
    }
    mysub = mysub[order(as.character(mysub$genes),decreasing=T),]
    circ <- circle_dat(mysub, genelist)
    pdf(paste(outdir,"BP.pdf",sep=""),width=20,height=10)
    GOBar(subset(circ, category == 'BP'),zsc.col = c('red','white','blue'))
    dev.off()
    pdf(paste(outdir,"MF.pdf",sep=""),width=20,height=10)
    GOBar(subset(circ, category == 'MF'),zsc.col = c('red','white','blue'))
    dev.off()
    pdf(paste(outdir,"CC.pdf",sep=""),width=20,height=10)
    GOBar(subset(circ, category == 'CC'),zsc.col = c('red','white','blue'))
    dev.off()
    #write XX_circ
    CC_circ = unique(subset(circ,category=="CC",select=c("category","ID","term","zscore","adj_pval")))
    BP_circ = unique(subset(circ,category=="BP",select=c("category","ID","term","zscore","adj_pval")))
    MF_circ = unique(subset(circ,category=="MF",select=c("category","ID","term","zscore","adj_pval")))
    write.table(MF_circ,file=paste(outdir,"MF_circ.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
    write.table(BP_circ,file=paste(outdir,"BP_circ.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
    write.table(CC_circ,file=paste(outdir,"CC_circ.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
    
    #Ramigo
    library(RamiGO)
    #set plot parameter
    pcolors <- c("white","yellow")
    psplit <-c(0.1,0.01,0.001,0.0001)
    zcolors = c("red","skyblue")
    zsplit <- c(8,4,2,1,0,-1,-2,-4)
    #plot
    if(length(BP_circ$ID)!=0){a = getAmigoTree(BP_circ$ID,pvalues=BP_circ$zscore,pcolors=zcolors,psplit=zsplit,filename=paste(outdir,"BPzscore",sep=""))}
    if(length(MF_circ$ID)!=0){a = getAmigoTree(MF_circ$ID,pvalues=MF_circ$zscore,pcolors=zcolors,psplit=zsplit,filename=paste(outdir,"MFzscore",sep=""))}
    if(length(CC_circ$ID)!=0){a = getAmigoTree(CC_circ$ID,pvalues=CC_circ$zscore,pcolors=zcolors,psplit=zsplit,filename=paste(outdir,"CCzscore",sep=""))}
    if(length(BP_circ$ID)!=0){a = getAmigoTree(BP_circ$ID,pvalues=BP_circ$adj_pval,pcolors=pcolors,psplit=psplit,filename=paste(outdir,"BPpval",sep=""))}
    if(length(MF_circ$ID)!=0){a = getAmigoTree(MF_circ$ID,pvalues=MF_circ$adj_pval,pcolors=pcolors,psplit=psplit,filename=paste(outdir,"MFpval",sep=""))}
    if(length(CC_circ$ID)!=0){a = getAmigoTree(CC_circ$ID,pvalues=CC_circ$adj_pval,pcolors=pcolors,psplit=psplit,filename=paste(outdir,"CCpval",sep=""))}
}

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

    
    outdir	outprefix (defalut:./goplot)
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
outdir <- ifelse(is.na(Args[3]),"goplot",Args[3])
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
#colnames(genelist)[5] <- "logFC"#GOplot format
universeData = read.table(dbstatsFile,header=F,col.names=c("go_id","evidence","gene_id"))
goData = universeData[universeData$gene_id %in% genelist$ID,]
params = GOstat_and_output()
# circle_dat function requires NOISeq-like genelist, hance,
# only on that circumstances, can we plot.
#plot
if(!ifonlylist)
{
    # The return value of prepare_for_plot,
    # and the input parameters of plot is under determined
    prepare_for_plot(params)
    plot()
}
