#!/usr/bin/env python

import sys, os, time, itertools
from multiprocessing import Pool
from collections import OrderedDict as Ordic

from RNAseqpipe.progsuit import Prog_Rsp, log, addends, get_filename
from RNAseqpipe.get_gene_length import len_for_Rsp

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
    """
    run verse and merge
    files: bam/sam
    outfiles: [countfile0,countfile1,]
    merged file write to mpath,
    each count file write to outpath
    """
    def _getpath(files,outpath,outfiles):
        mypath = []
        for i in range(len(outfiles)):
            myfile = files[i]
            try:
                path = outpath[i]
            except:
                path = os.path.dirname(myfile)
            mypath.append(path)
        return mypath

    assert isinstance(files,list)
    if outpath and len(files) != len(outpath):
        log.error("outpaths should be as many as bam/sam files are")
    outfiles = [None for i in range(len(files))]
    mypath = _getpath(files,outpath,outfiles)
    with Pool(p) as pool:
        outfiles = pool.starmap(verse,zip(itertools.repeat(conf),files,mypath,itertools.repeat(gff),itertools.repeat(silence)))
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

def main(argv):
    import argparse
    from RNAseqpipe.progsuit import Configuration
    from RNAseqpipe.run_RNAseqpipe import BASE_CONF

    parser = argparse.ArgumentParser(description="Count with verse")
    parser.add_argument('infile',nargs='+',help="input bam files")
    parser.add_argument('--gtf',required=True,help="gtf file")
    #parser.add_argument('-c','--conf',required=True,help="conf file")
    parser.add_argument('-c','--conf',help="conf file")
    parser.add_argument('-o','--outfile',required=True,help="output file")
    args = parser.parse_args(argv[1:])
    ## set to info level
    log.setLevel(20)
    myconf = Configuration(args.conf, base_conf=BASE_CONF)
    versepip(myconf,args.infile,args.outfile,args.gtf,silence=False)

if __name__ == '__main__':
    import sys
    main(sys.argv)
