#!/usr/bin/env python

'''
_version=0.0.1
'''
import sys, os, time, itertools
from multiprocessing import Pool
from collections import OrderedDict as Ordic

from RNAseqpipe.progsuit import Prog_Rsp, log, addends, get_filename
from RNAseqpipe.get_gene_length import len_for_Rsp

SALMON=os.path.dirname(sys.argv[0])+"/expression/salmon_quantify.sh"
def salmonpip(conf,myfq1,myfq2,fqnames,outdir,gff,p=8):
    #run salmon
    args = [SALMON,gff,conf["all"]["genome"],outdir,p]
    for i in range(len(myfq1)):
        args.extend((fqnames[i],myfq1[i],myfq2[i]))
    order = " ".join(map(str,args))
    #os.system(order)
    subprocess.call([order],shell=True)

def main(argv):
    import argparse
    from RNAseqpipe.progsuit import Configuration
    from RNAseqpipe.run_RNAseqpipe import BASE_CONF

    parser = argparse.ArgumentParser(description="Count with salmon")
    parser.add_argument('infile',nargs='?',help="Input RNAseqpipe data file.")
    parser.add_argument('--gtf',required=True,help="gtf file")
    #parser.add_argument('-c','--conf',required=True,help="conf file")
    parser.add_argument('-c','--conf',help="conf file")
    parser.add_argument('-o','--outfile',required=True,help="output file")
    args = parser.parse_args(argv[1:])
    ## set to info level
    log.setLevel(20)
    myconf = Configuration(args.conf, base_conf=BASE_CONF)
    # Parse group data
    # Group data are loaded in RNAseqpipe parser.
    mygroup_data = Group_data(args.group_data)
    fqnames, myfq1, myfq2 = mygroup_data.get2pip()

    salmonpip(myconf,args.infile,args.outfile,args.gtf,silence=False)
    #args=[SALMON]
    #args.extend(argv)
    #subprocess.call(map(str,args))
    
if __name__ == '__main__':

    import sys
    main(sys.argv[1:])
