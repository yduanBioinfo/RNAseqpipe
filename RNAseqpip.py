#!/usr/bin/env python
'''
_version = 0.2.3
wrapper for RNAseqpip programs
'''
from __future__ import print_function
import sys, re, os, copy
from threading import Thread
import subprocess
import logging

import hisat, cufflinks, NOISeq, htseq
import funcAnnot.b2gprog.GOannot as GOannot
import funcAnnot.kegg.KOannot as KOannot
from funcAnnot.kegg.annot2go_stat import annot2go_stat
from funcAnnot.annot_from_db import annot
from funcAnnot.kegg.koenrich.keggenrich import enrich as koenrich
from progsuit import Configuration, Group_data, getAbsPath, matchpath, log
from get_gene_length import len_for_Rsp
from get_geneids import get_geneids
    
logging.basicConfig()
FILEPATH=os.path.realpath(__file__)
FILEDIR=os.path.dirname(FILEPATH)
CPTDE="cptDE.py"
SEQ2EXP="seq2exp.py"
FUNC="func.py"

def run_subp(argv,program):
    #program shoule be in {CPTDE,SEQ2EXP,FUNC}
    del argv[0]
    argv[0] = FILEDIR+"/"+program
    subprocess.call(argv)

if __name__ == '__main__':

    import argparse, sys
    
    parser = argparse.ArgumentParser(description='RNA-seq analyse pip')
    parser.add_argument('program',help='all for whole pip/ali for alignment/ass for assembly/\
    cptDE for compute different expression',choices=['all','cptDE','seq2exp','func'])
    args=parser.parse_args(sys.argv[1:2])

    if args.program == 'cptDE':
        run_subp(sys.argv,CPTDE)        

    if args.program == 'seq2exp':
        run_subp(sys.argv,SEQ2EXP)        
        
    if args.program == 'func':
        run_subp(sys.argv,FUNC)        
        
