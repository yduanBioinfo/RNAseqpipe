#!/usr/bin/env python3

import sys, itertools, os
from multiprocessing import Pool
p = 16
script="/home/yduan/script/shell/qc_trim_rna_tmp.sh"

def split_fqs(input_fqs):
    """[fq1_1, fq1_2, fq2_1, fq2_2, fq3_1, fq3_2] -->
    [fq1_1,fq2_1, fq3_1], [fq1_2, fq2_2, fq3_2]"""
    fq1 = []
    fq2 = []
    for i in range(0, len(input_fqs), 2):
        fq1.append(input_fqs[i])
        fq2.append(input_fqs[i+1])
    return fq1, fq2

def run_qc(fq1, fq2, outdir='.'):
    order = " ".join(['bash', script, fq1, fq2])
    os.system(order)

def qc(input_fqs, outpath='.', maxp=20):
    #outfiles = []
    #def add_trans(path):
    #    return path+"/transcripts.gtf"

    #if not outpath:
    #    outpath = map(os.path.dirname,bamfs)
    #outfiles = map(add_trans,outpath)

    fq1, fq2 = split_fqs(input_fqs)
    pool = Pool(maxp)
    pool.starmap(run_qc, zip(fq1, fq2, itertools.repeat(outpath)))
    #pool.map(cufflink_star,zip(itertools.repeat(myconf),bamfs,outpath,itertools.repeat(silence)))
    #return outfiles

qc(sys.argv[1:],'./')

