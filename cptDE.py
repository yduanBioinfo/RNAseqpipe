#!/usr/bin/env python

from __future__ import print_function
import sys, os

import RNAseqpipe.NOISeq
from RNAseqpipe.progsuit import Configuration, Group_data, getAbsPath, matchpath
from RNAseqpipe.run_RNAseqpipe import add_arguments
from RNAseqpipe.get_gene_length import len_for_Rsp

def cptDE(myconf,ali_path,mygroup_data,expressionf,samples_table,template_gff):
    #args.expressionf args.samples_table args.template_gff
    
    if not expressionf or not samples_table:
        sys.stderr.write("expression file and samples_table file must be specific when \
running cptDE.\n")
        #raise StandardError
    expressionf = os.path.abspath(expressionf)
        
    #mergedfile : the gff file used in getting expression level(count or others)
    #mergedfile offers gene length for normalization
    if template_gff:
        mergedfile = template_gff
    else:
        mergedfile = myconf.get("all").get("gff")
    #mergedfile = open("/home/yduan/data/growth/ruibo/pip/version2/tmp/merged_asm/merged.gtf")
    try:
        mergedfpath = os.path.abspath(mergedfile)
        mergedfile = open(mergedfile)#len_for_Rsp required
    except:
        mergedfpath = os.path.abspath(mergedfile.name)#already a file object
    if not mergedfile:
        sys.stderr.write("gff file is needed for gene length\n")
        raise TypeError
    if myconf.get("gff") and myconf.get("gff").get("id_attr"):
        myt_attr = myconf.get("gff").get("id_attr")
    length = len_for_Rsp(mergedfile,t_attr=myt_attr)
    mergedfile.close()
    #mergedfile should be opened file
    mycounts = NOISeq.noi_counts(myconf,expressionf,mygroup_data.get_g_groups(),\
    mygroup_data.get_g_header(),sample_info=samples_table,qcreport=True,outpath=ali_path,\
    length=length)
    #mycounts :[n,rpkm,tmm,uqua]*[null,up,down]
    #mycounts :NOISeq inference for count data
    uqua_null,uqua_up,uqua_down = mycounts[3]
    print(mycounts)
    
    return uqua_null,uqua_up,uqua_down

def main(argv):

    import argparse, sys
    
    parser = argparse.ArgumentParser(description='cptDE program')
    parser.add_argument('-s','--samples_table',help='cuffnorm samples.table file'\
    ,nargs='?')
    parser.add_argument('-e','--expressionf',help='when chosing cptDE, this option is required.'\
    ,nargs='?')
    parser.add_argument('-t','--template_gff',help='when chosing cptDE or func,\
    this option is used.',nargs='?')
    add_arguments(parser)
    args = parser.parse_args(argv[1:])

    if len(argv) == 1:
        parser.print_usage()
        sys.exit(1)

    myconf = Configuration(args.conf)
    ali_path = getAbsPath(args.outpath)#home path for alignment results
    mygroup_data = Group_data(args.group_data)
        
    nullexp,upexp,downexp = cptDE(myconf,ali_path,mygroup_data,args.expressionf,args.samples_table,args.template_gff)

if __name__ == '__main__':

    import sys

    main(sys.argv)
