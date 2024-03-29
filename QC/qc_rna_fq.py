#!/usr/bin/env python3

import sys, os, copy, itertools, glob
import subprocess, pkg_resources
from subprocess import PIPE
from collections import OrderedDict as Ordic
from multiprocessing import Pool

from RNAseqpipe.progsuit import Configuration, Prog_Rsp, log

def split_fqs(input_fqs):
    """[fq1_1, fq1_2, fq2_1, fq2_2, fq3_1, fq3_2] -->
    [fq1_1,fq2_1, fq3_1], [fq1_2, fq2_2, fq3_2]"""
    fq1 = []
    fq2 = []
    for i in range(0, len(input_fqs), 2):
        fq1.append(input_fqs[i])
        fq2.append(input_fqs[i+1])
    return fq1, fq2

def fastqc_one_sample(myconf, fq, outdir):
    progname = 'fastqc'
    order1 = Ordic([('-o', outdir), (fq, '')])
    order2 = {}
    
    prog = Prog_Rsp(myconf, progname, order1, order2)
    prog.run()

def run_multiQC(myconf, fdir, sdir):
    '''sdir: summary directory,
    fdir: fastqc results directory.
    '''
    progname = 'multiqc'
    order1 = Ordic([('-o', sdir), (fdir, '')])
    order2 = {}
    
    prog = Prog_Rsp(myconf, progname, order1, order2)
    prog.run()

def fastqc_multqc(myconf, fqs, outdir, maxp=20): 
    """Run fastqc and multiQC"""
    fdir = os.path.join(outdir, 'fastqc')
    sdir = os.path.join(outdir, 'multiQC')
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(fdir, exist_ok=True)
    os.makedirs(sdir, exist_ok=True)
    pool = Pool(maxp)
    pool.starmap(fastqc_one_sample, zip(itertools.repeat(myconf), fqs, itertools.repeat(fdir)))
    # Run multiQC
    run_multiQC(myconf, fdir, sdir)

def trim_one_sample(myconf, fq1, fq2, outdir):
    # Set cores to 1
    cores = 1
    os.makedirs(outdir, exist_ok=True)
    # IlluQC
    filt_dir = os.path.join(outdir, 'filt_')
    os.makedirs(filt_dir, exist_ok=True)
    base1 = os.path.basename(fq1)
    base2 = os.path.basename(fq2)

    progname = 'IlluQC_PRLL.pl'
    order1 = Ordic([('-pe', fq1), (fq2, myconf['all'].get('adaptor')), 
        ('A', ''), ('-o', filt_dir)])
    order2 = Ordic([('-c', cores)])
    
    prog = Prog_Rsp(myconf, progname, order1, order2)
    prog.run()
    
    # AmbiguityFiltering
    qc1 = os.path.join(filt_dir, base1+'_filtered')
    qc2 = os.path.join(filt_dir, base2+'_filtered')
    progname = 'AmbiguityFiltering.pl'
    order1 = Ordic([('-i', qc1), ('-irev', qc2), ('-t3', '-t5')])
    order2 = {}
    
    prog = Prog_Rsp(myconf, progname, order1, order2)
    prog.run()

    # Mv results to clean_dir
    # And clean directory
    clean_dir = os.path.join(outdir, 'clean')
    fq1_clean = os.path.join(clean_dir, base1.split('.')[0]+'.clean_read1.fq')
    fq2_clean = os.path.join(clean_dir, base2.split('.')[0]+'.clean_read2.fq')
    os.makedirs(clean_dir, exist_ok=True)
    os.remove(qc1)
    os.remove(qc2)
    os.rename(qc1+'_trimmed', fq1_clean)
    os.rename(qc2+'_trimmed', fq2_clean)
    os.system('gzip {}'.format(fq1_clean))
    os.system('gzip {}'.format(fq2_clean))

def qc(myconf, input_fqs, outdir='.', maxp=20):
    # fastqc before control
    report_dir = os.path.join(outdir, 'report_before_qc')
    fastqc_multqc(myconf, input_fqs, report_dir, maxp)
    # Filt out bad sequences
    fq1, fq2 = split_fqs(input_fqs)
    pool = Pool(maxp)
    pool.starmap(trim_one_sample, zip(itertools.repeat(myconf), fq1, fq2, itertools.repeat(outdir)))
    # fastqc after control
    report_dir = os.path.join(outdir, 'report_after_qc')
    clean_fqs = glob.glob(os.path.join(outdir, 'clean','*fq.gz'))
    fastqc_multqc(myconf, clean_fqs, report_dir, maxp)

def main(argv):
    import argparse
    from RNAseqpipe.progsuit import Configuration
    from RNAseqpipe.run_RNAseqpipe import BASE_CONF

    parser = argparse.ArgumentParser(description="Quality control for illumina fastq files.")
    parser.add_argument('infile', nargs='+', help="input fastq files")
    parser.add_argument('-c', '--conf', help="Please set QC conf. Example: confs/QC.conf")
    parser.add_argument('-t', '--threading', type=int, help="threading used for computing")
    parser.add_argument('-o','--outdir',required=True,help="output directory")
    args = parser.parse_args(argv[1:])
    ## set to info level
    log.setLevel(20)

    qc_conf = args.conf if args.conf else \
             pkg_resources.resource_filename('RNAseqpipe', 'confs/QC.conf')
    myconf = Configuration(qc_conf, base_conf=BASE_CONF)
    qc(myconf, args.infile, args.outdir, args.threading)

if __name__ == '__main__':
    import sys
    main(sys.argv)
