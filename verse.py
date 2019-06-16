#!/usr/bin/env python

import sys, os, time, threading
from collections import OrderedDict as Ordic
from progsuit import Prog_Rsp, log, addends, get_filename
from get_gene_length import len_for_Rsp

mythreads = []#threading running
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

    def __init__(self,conf,file,path,gff,silence=False,tdindx=-1,\
    name='',outfiles=[]):

        super(Mythread,self).__init__()
        self.conf = conf
        self.file = file
        self.path = path
        self.gff = gff
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
            self.child(w)

    def child(self,w):#run verse

        fcount = verse(self.conf,self.file,self.path,self.gff,self.silence)
        #fcount: count file path (string)
        os.write(w,fcount)
        os._exit(0)

    def parent(self,r):

        info = ''
        data = os.read(r,64)
        while data:
            info += data
            data = os.read(r,64)

        write_ed(self.tdindx,info,self.name,self.outfiles)

def verse(conf,bam,outpath,gff,silence=False):
    
    #run one verse program

    progname = "verse"
    filename = get_filename(bam)
    outfile = addends(outpath)+filename
    order1 = {"-o":outfile,"-a":gff,bam:""}
    order2 = {}

    prog = Prog_Rsp(conf,progname,order1,order2,silence)
    prog.run()
    return outfile+".exon.txt"

#merge count file
# lfile: length file
def catcount(conf,files,outfile,lfile=None,silence=False):

    progname = "count_merge"
    order1 = Ordic([(file,"") for file in files])
    order2 = Ordic([("-o",outfile),("-f","")])
    if lfile:
        order2["-l"]=lfile
    prog = Prog_Rsp(conf,progname,order1,order2,silence)
    return prog.run()

def get_length(conf,outfile,gff,silence=False):
    #merge count file

    progname = "get_length"
    order1 = Ordic([(gff,"")])
    order2 = Ordic([("-o",outfile)])
    prog = Prog_Rsp(conf,progname,order1,order2,silence)
    return prog.run()

def versepip(conf,files,mpath,gff,outpath=None,p=8,silence=False):

    #run verse and merge
    #files: bam/sam
    #outfiles: [countfile0,countfile1,]
    #merged file write to mpath,
    #each count file write to outpath

    assert isinstance(files,list)
    if outpath and len(files) != len(outpath):
        log.error("outpaths should be as many as bam/sam files are")
    outfiles = [None for i in range(len(files))]

    for i in range(len(outfiles)):
        myfile = files[i]
        try:
            path = outpath[i]
        except:
            path = os.path.dirname(myfile)

        tname = "thread%d"%i
        mythread = Mythread(conf,myfile,path,gff,silence,i,name=tname,outfiles=outfiles)
        mythreads.append(tname)
        mythread.start()
        waittd(p)
    #waiting for all threads done
    waittd(0)

    outfile = addends(mpath)+time.strftime("merged%y_%m_%d_%H.count",time.localtime())
    # Add length to outfile
    
    # Block get length... the core script is in alter need. 
    # The error info
    #File "/home/yduan/script/python/bio/seq/base.py", line 235, in __iterGff
    #yield Gff_rec(line.split(self.sep),self.fm)
    #File "/home/yduan/script/python/bio/seq/base.py", line 290, in __init__
    #self.format_attr(fm)
    #File "/home/yduan/script/python/bio/seq/base.py", line 321, in format_attr
    #self.format_normal()
    #File "/home/yduan/script/python/bio/seq/base.py", line 311, in format_normal
    #self.gene_id = self.attr["gene_id"]
    #KeyError: 'gene_id'

    #lenfile = addends(mpath)+time.strftime("merged%y_%m_%d_%H.length",time.localtime())
    #get_length(conf,lenfile,gff)
    lenfile = None
    catcount(conf,outfiles,outfile,lfile=lenfile)

    return outfile
