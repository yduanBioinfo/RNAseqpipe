#!/usr/bin/env python3

import os
import subprocess
from glob import glob
from bio.seq.base import Fasta, write_fasta_o

'''
_version=0.0.2
'''

def make_path(outpath):
    if outpath:
        outpath = os.path.abspath(outpath)
        if not outpath.endswith("/"):
            outpath = outpath+"/"
            try:
                os.mkdir(outpath)
            except OSError:
                pass
    else :
        outpath = ""
    return outpath

def subset_fa(infa,ids,outfa):
    with open(outfa,'w') as outfa:
        for rec in Fasta(infa,if_trim=True):
            if rec.name in ids:
                write_fasta_o(outfa,rec)

def get_filename(path):
    """Get basename with out suffix."""
    return os.path.splitext(os.path.basename(path))[0]
    
def get_fa_IDs(infile,trim_id=True):
    """ If trim_id, only keeps the first field. eg, the name of ">aaa bbb ccc" is aaa.
    """
    ids = []
    for seq in Fasta(infile,if_trim = trim_id):
        ids.append(seq.name)
    return frozenset(ids)

def get_IDs(infile,header=False,col=0,writefile=False,outpath=""):

    '''infile : NOISeq infile [defualt]
       read one column into list.
    '''
    
    if isinstance(infile,str):
        fp = open(infile)
    else:
        fp = infile
    IDs = []
    if header:
        fp.readline()
    for each in fp:
        tmp = each.split("\t")
        if tmp and tmp[col].strip() not in IDs:
            IDs.append(tmp[col].strip())
    if writefile:
        IDfile = open(outpath+get_filename(file)+".ID",'w')
        IDfile.write("\n".join(IDs))
        IDfile.close()
        
    return frozenset(IDs)

def get_seq(gtf_file,genome,outfile):
    subprocess.call(["gffread",gtf_file,"-g",genome,"-w",outfile])
    
def blastx(genome,fasta,outfile,blastx="/home/yduan/soft/bio/seq/ncbi-blast-2.2.27+/bin/blastx"):
    subprocess.call([blastx,"-query",fasta,"-db",genome,"-out",outfile,"-outfmt","5","-evalue","0.000001","-max_target_seqs","5","-num_threads",str(thread),"-show_gis"])
    
def b2g(xmlfile,outfile):
    subprocess.call(["b2g","-in",xmlfile,"-annot","-out",outfile])

def fasta2annot(fastafile,db,outpath="",thread=8):
    '''
    fasta -> xml -> annot
    xml name:fastafile+.xml
    annot file name: fastafile+.annot
    '''
    
    xmlfile = outpath+get_filename(fastafile)+".xml"
    blastx(db,fastafile,xmlfile,thread)
    annotfile = xmlfile[:-4]#del ".xml" 
    b2g(xmlfile,annotfile)#blast2go out annotfile+".annot"
    return annotfile+".annot"

def annot_gff(ingff,genome,dbs,outpath,thread):
    outfile=outpath+get_filename(ingff)+".fa"
    get_seq(ingff,genome,outfile)
    annot_fa(outfile,dbs,outpath,thread)
    subprocess.call(["rm",outfile])

def clean_tmp(outpath,fa=True,annot=True,xml=False):
    clean_lst = ["rm"]
    if fa:
        clean_lst.extend(glob(outpath+"tmp*.fa"))
    if annot:
        clean_lst.extend(glob(outpath+"tmp*.annot"))
    if xml:
        clean_lst.extend(glob(outpath+"tmp*.xml"))
    if len(clean_lst) > 1:
        subprocess.call(clean_lst)

def annot_fa(infa,dbs,outpath,thread):
    '''
    fasta -> xml -> annot -> IDs -> leftIDs -> unannoted fasta (-> xml)
    fasta name: tmp+index_of_db+.fa
    xml name: tmp+index_of_db+.xml
    annot file name: tmp+index_of_db+.annot
    Result file: final.annot
    '''
    end_annot = outpath+"final.annot"
    total_id = get_fa_IDs(infa)
    for i in range(len(dbs)):
        db = dbs[i]
        tmp_fa = outpath+"tmp"+str(i)+".fa"
        subset_fa(infa,total_id,tmp_fa)
        annote = fasta2annot(tmp_fa,db,outpath,thread)
        annote_id = get_IDs(annote)
        total_id = total_id - annote_id
        infa = tmp_fa
    outfiles = glob(outpath+"tmp*.annot")
    args = ["cat"]
    args.extend(outfiles)
    args.append(">"+end_annot)
    subprocess.call(" ".join(args),shell=True)
    clean_tmp(outpath)

def main(argv):
    
    import argparse, sys, os
    
    parser = argparse.ArgumentParser(description='blast2go pip')
    parser.add_argument('--infile',help='Fasta or gff file. Genome is in need, when gff is provide',nargs='?')
    parser.add_argument('-s','--genome',nargs='?',help='if genome fasta file is provide, infile should be gff file')
    parser.add_argument('-d','--dbs',nargs='+',help='protein dbs, .faa file')
    parser.add_argument('-t','--thread',nargs='?',help='Number of threading for blastx',default=8)
    parser.add_argument('-o','--outpath',nargs='?',help='out file path',default='./b2g')
    args = parser.parse_args(argv[1:])    
    
    outpath = make_path(args.outpath)
    if args.genome:
        annot_gff(args.infile,args.genome,args.dbs,outpath,args.thread)
    else:
        annot_fa(args.infile,args.dbs,outpath,args.thread)
    
if __name__ == '__main__':

    import sys
    
    main(sys.argv)
    
