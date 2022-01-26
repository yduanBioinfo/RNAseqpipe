#!/usr/bin/env python

#_version=0.0.2

import sys, re 
import numpy as np
from collections import OrderedDict as Ordic
from RNAseqpipe.progsuit import Prog_Rsp, log
try:
    from pyper import *
except:
    log.warning("pyper has not installed. Install this package before running NOISeq analysis.")

def is_number(mystr):

    return re.match(r'\d+\.?\d*$',mystr)#16/8/18 add '$'

def startswithd(mystr):

    return re.match(r'\d+.*',mystr)
    
def r_str(mystr):
    #convert all parameters to str, str to "str".
        
    if mystr in ["FALSE","TRUE","F","T","NULL"]:
        return mystr
    try:
        if is_number(mystr):
            return mystr
    except:#mystr is not a type of string
        return str(mystr)
    return '\"'+mystr+'\"'
    
def writerr(comm,verbose=True):
    
    if verbose:
        sys.stderr.write(comm)
    
def make_group_dic(sample_info,has_header=True,sep="\t"):
    #outdic {q_id: sample_name}
    #cufflinks required
    open_file = False
    if isinstance(sample_info,str):
        sample_info = open(sample_info)
        open_file = True
        
    outdic = Ordic()
    if has_header:
        header = sample_info.readline()
    for eachline in sample_info:
        tmp = eachline.split(sep)
        name = tmp[1].split("/")[-2]
        outdic[tmp[0]] = name
        
    if open_file:
        sample_info.close()
        
    return outdic
        
def sorted_group(groups,group_dic,indx=0,has_header=False,sep="\t"):
    #group_dic {q_id:name}(q_id: cuffnorm out sample id;name: sample name)
    #cufflinks required
    #indx: index of name column
    
    if has_header:
        mygroups = groups[1:]
    else:
        mygroups = groups[:]
    names = map(lambda x:x[indx],mygroups)
    out_groups = []
    for key, value in group_dic.items():
        out_groups.append(mygroups[names.index(value)])
    return out_groups

def add_X(mystr,p="X"):
    #add "X" as prefix for R when comes across a digit or a string start with digit
    #return a string
    #p:prefix

    try:
        if startswithd(mystr):
            return p+mystr
    except:#mystr is not a type of string
        sys.stderr.write("Warning:add_X using in an unexpect way!!!")
        return p+str(mystr)
    return mystr

def altgroup(mygroup,indx=0):
    
    return map(lambda x:x[0:indx]+(add_X(x[indx]),)+x[indx+1:],mygroup)
    
def NOISeq(data,groups):
    pass
    
def NOISeqBIO(data,groups,header,vb,length,k=0.01, norm = "rpkm",lc=1, pr = 20, adj = 1.5, \
plot = 'FALSE', a0per = 0.9, random_seed = 12345,filter = 1,cpm = 0.05,\
outname = "NOIbio",factor="",q=0.8,sample_info=None,qcreport=False,qcfile='NULL',outpath=''):
    #header: header of groups
    #r in R changed to pr/random.seed changed to random_seed
    #sample_info: cuffnorm out file: samples.table[cufflinks parameter]
    #qcfile can only write to work path, due to QCreport function limits.!!!
    #vb : verbose
    print("start!!!")
    sys.stdout.flush()
    if outpath and not outpath.endswith("/"):
        outpath = outpath+"/"
    sep = "\t"
    outname = outpath+outname
        
    def make_header(header,samples='samples',strc='S18'):
        #strc: structure type of array in np
        outheader = header[:]
        outheader.insert(0,samples)
        for i in range(len(outheader)):
            outheader[i] = (outheader[i],strc)
        return outheader
    
    def r_str(mystr):
        #convert all parameters to str, str to "str".
        
        if mystr in ["FALSE","TRUE","F","T","NULL"]:
            return mystr
        try:
            if is_number(mystr):
                return mystr
        except:#mystr is not a type of string
            return str(mystr)
        return '\"'+mystr+'\"'

    #sort group
    #print(is_number(lc))
    group_dic = None
    
    if sample_info:
        group_dic = make_group_dic(sample_info)        
    if group_dic:
        groups = sorted_group(groups,group_dic)
    #1
    print("groups:")
    print(groups)
    groups = altgroup(groups)
    #x
    print("groups:")
    print(groups)
    sys.stdout.flush()
    if factor not in header:
        factor = header[-1]
        
    data, k, norm, lc, pr, adj, plot, a0per, random_seed, filter, cpm, factor, q, sep= \
    map(r_str,[data, k,norm, lc, pr, adj, plot, a0per, random_seed, filter, cpm, factor, q, sep])
    
    r = R()
    writerr(r("library(NOISeq)"),vb)
    writerr(r("mytable<-read.table(%s,header=T)"%data),vb)
    writerr(r('row.names(mytable)<-mytable[,1]\n\
    mytable<-mytable[,-1]'),vb)
    r.myfactors = np.array(groups,make_header(header))
    #2
    print("myfactors:")
    print(r("myfactors"))
    sys.stdout.flush()
    r.mylength = "NULL"
    if length:
        r.length = np.array(length,[("gene",'|S64'),("length",'i4')])
        #r.length = np.array(length,[("gene",'U'),("length",'i4')])
        writerr(r('mylength<-as.integer(as.vector(length[,2]))\n\
        names(mylength)<-length[,1]'),vb)
        writerr(r('head(mylength)'),vb)
        writerr(r('mylength<-subset(mylength,names(mylength) %in% rownames(mytable))'),vb)
        writerr(r('head(mylength)'),vb)
    #3
    print("myfactors[,1]:")
    print(r("myfactors[,1]"))
    
    print("old mytable:")
    print(r("head(mytable)"))
    sys.stdout.flush()
    writerr(r('mytable<-mytable[as.vector(myfactors[,1])]'),vb)
    #4
    print("mytable:")
    print(r("head(mytable)"))
    sys.stdout.flush()
    if norm != '"n"':
        print("norm:%s" % norm)
        if norm == '"tmm"':
            writerr(r('mynorm<-%s(mytable,long=mylength,lc=%s,k=%s)'%("tmm",lc,k)),vb)
            #writerr(r('mynorm<-%s(mytable)'%"tmm"),vb)
        if norm == '"rpkm"':
            writerr(r('mynorm<-%s(mytable,long=mylength,lc=%s,k=%s)'%("rpkm",lc,k)),vb)
        if norm == '"uqua"':
            writerr(r('mynorm<-%s(mytable,long=mylength,lc=%s,k=%s)'%("uqua",lc,k)),vb)
        print(r('head(mynorm)'))
        #writerr(r('mynorm<-%s(mytable,long=mylength,lc=%s,k=%s)'%(norm.strip('"'),lc,k)),vb)
        writerr(r('write.table(mynorm,file=%s,sep=%s,quote=FALSE)'%(r_str(outname+norm.strip('"')+".txt"),sep)),vb)
    writerr(r('mydata<-readData(data=mytable,factors=myfactors,length=mylength)'),vb)
    writerr(r('mynoiseqbio<-noiseqbio(mydata,k=%s,norm=%s,lc=%s,r=%s,adj=%s,plot=%s,\
    a0per=%s,random.seed=%s,filter=%s,cpm=%s,factor=%s)'%\
    (k,norm,lc,pr,adj,plot,a0per,random_seed,filter,cpm,factor)),vb)
    writerr(r('mynoiseq.deg = degenes(mynoiseqbio, q = %s, M = NULL)'%q),vb)
    writerr(r('write.table(mynoiseq.deg,file=%s,sep=%s,quote=FALSE)'%(r_str(outname+"null.txt"),sep)),vb)
    writerr(r('mynoiseq.deg = degenes(mynoiseqbio, q = %s, M = \"up\")'%q),vb)
    writerr(r('write.table(mynoiseq.deg,file=%s,sep=%s,quote=FALSE)'%(r_str(outname+"up.txt"),sep)),vb)
    writerr(r('mynoiseq.deg = degenes(mynoiseqbio, q = %s, M = \"down\")'%q),vb)
    writerr(r('write.table(mynoiseq.deg,file=%s,sep=%s,quote=FALSE)'%(r_str(outname+"down.txt"),sep)),vb)
    if qcreport:
        QCreport(r,'mydata',file=qcfile,factor=factor,path=outpath)
        QCreport(r,'mydata',file=qcfile,factor=factor,norm='TRUE',path=outpath)
    del r
    return outname+"null.txt", outname+"up.txt", outname+"down.txt"
    
def NOISeqBIO_count():

    pass
    
def noisb_fpkm(myconf,data,groups,header,verbose=True,length=None,k="NULL", norm = "n",lc=0, pr = 20,\
adj = 1.5, plot = 'FALSE', a0per = 0.9, random_seed = 12345,filter = 1,cpm = 1,\
outname = "NOIbio",factor="",q=0.95,sample_info='None',qcreport='False',qcfile='NULL',outpath=''):
    #NOISeqBIO_fpkm
    #length: length information of genes.(id,length) iterator is required
    
    order2 = {'k':k,'norm':norm,'lc':lc,'pr':pr,'adj':adj,'plot':plot,'a0per':a0per,\
    'random_seed':random_seed,'filter':filter,'cpm':cpm,'outname':outname,'factor':factor,\
    'q':q,'sample_info':sample_info,'qcreport':qcreport,'qcfile':qcfile,'outpath':outpath}
    mynoiseq = NOISeqbio(myconf,"noiseq",data,groups,header,verbose,length,order2=order2)#order2=order2
    return mynoiseq.run() 
    #return NOISeqBIO(data,groups,header,k,norm,lc,pr,adj,plot,a0per,random_seed,filter,cpm,\
    #outname,factor,q,sample_info,qcreport,qcfile,outpath)

def noisb_count(myconf,data,groups,header,verbose=True,length=None,k=0.5, norm = 'rpkm',lc=1, pr = 20,\
adj = 1.5, plot = 'FALSE', a0per = 0.9, random_seed = 12345,filter = 1,cpm = 1,\
outname = "NOIbio",factor="",q=0.95,sample_info='None',qcreport=True,qcfile='NULL',outpath=''):
    #NOISeqBIO_fpkm
    #length: length information of genes.(id,length) iterator is required
    order2 = {'k':k,'norm':norm,'lc':lc,'pr':pr,'adj':adj,'plot':plot,'a0per':a0per,\
    'random_seed':random_seed,'filter':filter,'cpm':cpm,'outname':outname,'factor':factor,\
    'q':q,'sample_info':sample_info,'qcreport':qcreport,'qcfile':qcfile,'outpath':outpath}
    mynoiseq = NOISeqbio(myconf,"noiseq",data,groups,header,verbose,length,order2=order2)#order2=order2
    result = mynoiseq.run() 
    #del mynoiseq
    return result
    
def noi_counts(myconf,data,groups,header,verbose=True,length=None,k=0.5,norm = 'rpkm',\
lc=1, pr = 20,adj = 1.5, plot = 'FALSE', a0per = 0.9,random_seed = 12345,filter = 1,\
cpm = 1,outname = "NOIbio",factor="",q=0.95,sample_info='None',qcreport=True,qcfile='NULL',outpath=''):
    #run n/fpkm/tmm/uqau normalization
    #NOISeq in configuration.txt should be empty
    
    def run():
        return noisb_count(myconf,data,groups,header,verbose,length,k,norm,lc,pr,\
        adj, plot, a0per, random_seed ,filter ,cpm ,outname ,factor ,q ,sample_info ,\
        qcreport, qcfile ,outpath)
    
    outfiles = []#save noisb_count result
    
    norm='n'
    lc=1
    qcreport = True
    outname = "n_"+outname
    outfiles.append(run())
    
    norm='rpkm'
    qcreport = False
    outname = "rpkm_"+outname.lstrip("n_")
    outfiles.append(run())
    
    norm='tmm'
    outname = "tmm_"+outname.lstrip("rpkm_")
    outfiles.append(run())
    
    norm='uqua'
    outname = 'uqua_'+outname.lstrip("tmm_")
    outfiles.append(run())
    
    return outfiles
    
def QCreport(r,input, file = 'NULL', samples = 'NULL', factor = 'NULL',norm = 'FALSE',\
verbose=True,path=''):
    #r: R obj
    import time
    
    file = time.strftime("QCreport%y_%m_%d_%H_%M_%S.pdf",time.localtime())
    if path:
        #path should endwith "/"
        file = path+file
    print(r('QCreport(%s,file=%s,samples=%s,factor=%s,norm=%s)'%(input,r_str(file),samples,factor,norm)))
        
class NOISeqbio(Prog_Rsp):

    def __init__(self,Conf,progname,data,groups,header,verbose=True,length=None,order1={},order2={},silence=False):
    
        self.ex_data = data
        self.groups = groups
        self.header = header
        progname = "NOISeq"
        super(NOISeqbio,self).__init__(Conf,progname,order1,order2,silence)
        self.vb = verbose
        self.length = length

    def run_order(self,order,name="unknown",silence=False):
    
        return eval(order)

    def gen_order(self,program,orders):
        #generate order for run 
        #orders is a dict

        def r_str(mystr):
            try:
                if is_number(mystr):
                    #number in ""
                    return eval(mystr)
            except:
                #exactlly number, int or float
                if isinstance(mystr,str):
                    raise TypeError
                return mystr
            if mystr in ["False","True","None"]:
                return mystr
            return "\""+mystr+"\""
                
        tmporder = ""
        for key, value in orders.items():
            tmporder += ","+str(key)+"="+str(r_str(value))
        outorder = "NOISeqBIO(self.ex_data,self.groups,self.header,self.vb,self.length"+tmporder+")"

        return outorder
        
def main(args):

    from RNAseqpipe.progsuit import Configuration, Group_data
    
    Conf = Configuration(open(args[0]))
    group_data = Group_data(open(args[1]))
    data = args[2]
    groups = group_data.get_g_groups()
    header = group_data.get_g_header()
    
    return noisb_fpkm(Conf,data,groups,header,qcreport=True)
    
if __name__ == '__main__':
    
    import sys
    
    print(main(sys.argv[1:]))
