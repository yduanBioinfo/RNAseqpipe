#!/usr/bin/env python

'''
'''

import sys, os, copy
import itertools
from collections import OrderedDict as Ordic
from progsuit import Configuration, Prog_Rsp
from multiprocessing import Pool

def stringtie(conf,bamf,outfile,silence=False,thre=0):
    #conf if configuration obj.
    #bamf: bam file [path]
    #thre : -p value
    progname = "stringtie"
    order1 = {bamf:"","-o":outfile}
    order2 = {}
    prog = Prog_Rsp(conf,progname,order1,order2,silence)
    return prog.run()

def stringtie_star(a_b):
    # convert f([1,2]) to f(1,2)
    return stringtie(*a_b)

def stringties(conf,bamfs,outpath=None,silence=False,maxp=20):
    #conf is a Configuration obj
    #bamfs: bam/sam files [list]
    #if outpath there isn't a specific outpath,
    #out transcripts.gtf should write in the same dir of bam file.
    myconf = conf
    outfiles = []
    def add_trans(path):
        return path+"/transcripts.gtf"

    if outpath and len(bamfs) != len(outpath):
        error_handle(1,"stringties")
    if not outpath:
        outpath = map(os.path.dirname,bamfs)
    outfiles = map(add_trans,outpath)

    pool = Pool(maxp)
    pool.map(stringtie_star,itertools.izip(itertools.repeat(myconf),bamfs,outfiles,itertools.repeat(silence)))
    return outfiles

def merge(conf,gtfs,outfile,silence=False):
    
    order1 = Ordic([("--merge",""),("-o",outfile)])
    order2 = Ordic()
    if conf["all"]["gff"]:
        order2["-G"] = conf["all"]["gff"]
    order2[" ".join(gtfs)] = ""
    merge = Prog_Rsp(conf,"stringtie",order1,order2,silence)
    merge.run()
    del merge
    return outfile
    
def error_handle(code,name="unknown"):

    if code == 1:
        sys.stderr.write("Error: outpaths should be as many as bam/sam files are")
    if code:
        sys.stderr.write("Error in %s step.\n" % name)
        sys.exit(20)

if __name__ == '__main__':            
   pass 
