#!/usr/bin/env python

'''
_version=0.0.1
'''
import sys, os, time, threading
from collections import OrderedDict as Ordic

from RNAseqpipe.progsuit import Configuration, Prog_Rsp

mythreads = []
threadLock = threading.Lock()

def write_ed(i,info,threadname,outfiles):
    #write end info of multiprocess result
    #i the ith thread 
        
    threadLock.acquire()
    outfiles[i] = info
    mythreads.remove(threadname)
    threadLock.release()    
        
def waittd(max_p=8):#waiting for all threads done
    
    while len(mythreads):
        if len(mythreads) < max_p:break
        time.sleep(1)
            
class Mythread(threading.Thread):

    def __init__(self,conf,file,path,gff,plotqa=True,silence=False,tdindx=-1,\
    name='',outfiles=[]):
    
        super(Mythread,self).__init__()
        self.conf = conf
        self.file = file
        self.path = path
        self.gff = gff
        self.qa = plotqa
        self.silence = silence
        self.tdindx = tdindx #index of thread
        self.name = name
        self.outfiles = outfiles
        
    def run(self):
    
        r,w = os.pipe()
        pid = os.fork()
        if pid:
            os.close(w)
            self.parent(r)
        else:
            os.close(r)
            self.child(w,self.qa)

    def child(self,w,qa):#run htseq-count and htseq-qa
        
        if qa:
            fqa = htseq_qa(self.conf,self.file,self.path,self.silence)
        fcount = htseq_count(self.conf,self.file,self.path,self.gff,self.silence)
        #fqa, fcount: qa,count file path (string)
        os.write(w,fcount)
        os._exit(0)

    def parent(self,r):

        info = ''
        data = os.read(r,64)
        while data:
            info += data
            data = os.read(r,64)
            
        write_ed(self.tdindx,info,self.name,self.outfiles)
    
def get_filename(path):
    
    return os.path.splitext(os.path.split(path)[1])[0]
    
def addends(value,ends="/"):

    if value.endswith(ends):
        return value
    else:
        return value+"/"
        
def htseq_count(conf,bam,outpath,gff,silence=False,):

    progname = "htseq_count"
    filename = get_filename(bam)
    order1 = Ordic([(bam,"")])
    order1[gff] = ""
    outfile = addends(outpath)+filename+".count"
    order2 = Ordic([(">",outfile)])
    
    prog = Prog_Rsp(conf,progname,order1,order2,silence)
    prog.run()
    return outfile
    
def htseq_qa(conf,bam,outpath,silence=False):

    progname = "htseq_qa"
    filename = get_filename(bam)
    outfile = addends(outpath)+filename
    order1 = {bam:"","-o":outfile}
    order2 = {}
    htseq_qa = Prog_Rsp(conf,progname,order1,order2,silence)
    htseq_qa.run()
    return outfile

def catcount(conf,files,outfile,silence=False):
    #merge count file 

    progname = "count_merge"
    order1 = Ordic([(file,"") for file in files])
    order2 = Ordic([("-o",outfile)])
    prog = Prog_Rsp(conf,progname,order1,order2,silence)
    return prog.run()
    
def htseqpip(conf,files,mpath,gff,outpath=None,plotqa=True,p=8,silence=False):
    #run htseq_count and merge
    #files: bam/sam
    #outfiles: [countfile0,countfile1,] don't put out qa file
    #for no more use of it.
    #merged file write to mpath,
    #each count file write to outpath
    
    assert isinstance(files,list)
    if outpath and len(files) != len(outpath):
        error_handle(1,"htseq-count")
        
    outfiles = [None for i in range(len(files))]
    
    for i in range(len(files)):
        file = files[i]
        try:
            path = outpath[i]
        except:
            path = os.path.dirname(file)
        
        tname = "thread%d"%i
        mythread = Mythread(conf,file,path,gff,plotqa,silence,i,name=tname,outfiles=outfiles)
        mythreads.append(tname)
        mythread.start()
        waittd(p)
    #waiting for all threads done
    waittd(0)
    
    outfile = addends(mpath)+time.strftime("merged%y_%m_%d_%H.count",time.localtime())
    catcount(conf,outfiles,outfile)
    #outfiles are single count file
    #outfile is the merged file which delete the lines startswith __
    
    return outfile
    
class HTcount(Prog_Rsp):

    def __init__(self,Conf,progname,order1={},order2={},silence=False,appendorder=""):
    
        super(HTcount,self).__init__(Conf,progname,order1,order2,silence)
        self.apdord = appendorder
        
    def run(self):
        #for outside use
    
        order = self.gen_order(self.progname,self.data)
        order += self.apdord
        sys.stderr.write(order+"\n")
        return self.run_order(order,self.progname,self.silence)

def error_handle(code,name="unknown"):

    if code == 1:
        sys.stderr.write("Error: outpaths should be as many as bam/sam files are")
    if code:
        sys.stderr.write("Error in %s step.\n" % name)
        sys.exit(20)

def test(argv):

    with open(argv[0]) as conf:
        myconf = Configuration(conf)
    bams = argv[1:4]#3 test file
    outpath = argv[4]
    gff = argv[5]
    print(htseqpip(myconf,bams,outpath,gff))

def main(argv):

    import argparse, sys, os
    
    parser = argparse.ArgumentParser(description='htseq pip(count each bam/sam and merge result)')
    parser.add_argument('infile',help='bam/sam files',nargs='+')
    parser.add_argument('-c','--conf',type=argparse.FileType('r'),nargs='?',help='configuration file')
    parser.add_argument('-g','--gff_file',nargs='?',help='gff file')
    parser.add_argument('-qa','--plotqa',default=False,action='store_true',\
    help='whether plot qa reports')
    parser.add_argument('-p','--proc',nargs='?',help='processer numbers',type=int,default=8)
    parser.add_argument('-o','--outpath',nargs='?',help='out file path',default='./')
    args = parser.parse_args(argv)    
    
    myconf = Configuration(args.conf)
    htseqpip(myconf,args.infile,args.outpath,args.gff_file,plotqa=args.plotqa,p=args.proc)
    
if __name__ == '__main__':

    import sys
    #htseq bamf gff outpath configuration
    #with open(sys.argv[4]) as conf:
    #    myconf = Configuration(conf)
    #htseq_count(myconf,sys.argv[1],sys.argv[3],gff=sys.argv[2])
    #htseq-qa bam outpath
    #htseq_qa(myconf,sys.argv[1],sys.argv[2])
    main(sys.argv[1:])
