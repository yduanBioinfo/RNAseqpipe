#!/usr/bin/env python

from __future__ import print_function
import sys, os
from threading import Thread

import RNAseqpipe.funcAnnot.b2gprog.GOannot as GOannot
import RNAseqpipe.funcAnnot.kegg.KOannot as KOannot
from RNAseqpipe.funcAnnot.kegg.annot2go_stat import annot2go_stat
from RNAseqpipe.funcAnnot.annot_from_db import annot
from RNAseqpipe.funcAnnot.kegg.koenrich.keggenrich import enrich as koenrich
from RNAseqpipe.progsuit import Configuration, Group_data, getAbsPath, matchpath
from RNAseqpipe.run_RNAseqpipe import log, add_arguments
from RNAseqpipe.get_gene_length import len_for_Rsp
from RNAseqpipe.get_geneids import get_geneids

def func_annot(myconf,DEfile,template_gff,outpath=None,godb=None,kodb=None):#outpath should be bug
    #functiong annot,GO and KEGG
    #enrichment
    #visualisation
    #DEfile: difference expression file
    ogodb = godb
    okodb = kodb#for soft link 16/8/17
    godb = godb if os.path.isfile(godb) else None
    kodb = kodb if os.path.isfile(kodb) else None
    
    if not outpath:
        outpath = getAbsPath("./func")
        
    def go_annot(godb):
    
        godir = getAbsPath(outpath+"/goannot")
        altdb = godir+"/mygodb.annot"
        annot2go_stat(godb,altdb)#add "ISA"
        print("altdb;DEfile: %s\t%s"%(altdb,DEfile))
        os.system("%s/funcAnnot/b2gprog/gostat_goplot_RamiGO_in.R %s %s %s"%\
        (os.path.dirname(os.path.realpath(__file__)),altdb,DEfile,godir))#outfile??
    
    def ko_annot(kodb):
        
        kodir = getAbsPath(outpath+"/koannot")
        koout = kodir+"/ko.annot"
        annot(DEfile,kodb,koout)
        kgmap = myconf.get("all").get("kgmap")
        print("enrich start")
        koenrich(koout,kodb,kgmap)
        print("enrich done")
        
    def annot_godb():#annot db and annot DE genes
        godb,godb1 = GOannot.blast2go(myconf,geneidf,godir,gff)#why none in last version?
        try:os.symlink(godb,ogodb)#16/8/17
        except:pass
        go_annot(godb)
        
    def annot_kodb():#annot db and annot DE genes
        kodb = KOannot.koannot(myconf,geneidf,kodir,gff)
        try:os.symlink(kodb,okodb)#16/8/17
        except:pass
        ko_annot(kodb)
        
    if not godb or not kodb:
        if template_gff:
            gff = template_gff
        else:
            gff = myconf.get("all").get("gff")
        gff = os.path.abspath(gff)
        if not gff:
            sys.stderr.write("gff file is needed for go/ko annot\n")
            raise TypeError
        gffobj = open(gff)#len_for_Rsp required

        geneidf = gff+"geneids"

        get_geneids(gff,geneidf)
        gffdir = os.path.dirname(gff)
        godir = gffdir+"/goannot"
        kodir = gffdir+"/koannot"

    if not godb:#go annot for whole transcriptom
        #gff are not annoted
        gotd = Thread(target=annot_godb)
    else:
        gotd = Thread(target=go_annot,args=(godb,))
    gotd.start()
        
    if not kodb:#ko annot
        kotd = Thread(target=annot_kodb)
    else:
        kotd = Thread(target=ko_annot,args=(kodb,))
    kotd.start()
        
    try:#waiting for db annot complete
        gotd.join()
        kotd.join()
    except:
        pass

def main(argv):    
    
    import argparse, sys
    
    parser = argparse.ArgumentParser(description='RNASeqpipe function annot program')
    parser.add_argument('-d','--DEfile',help='DEfile(NOISeq out). when chosing func,\
    this option is required.',nargs='?')
    parser.add_argument('-t','--template_gff',help='when chosing cptDE or func,\
    this option is used.',nargs='?')
    add_arguments(parser)
    args = parser.parse_args(argv[1:])

    if len(argv) == 1:
        parser.print_usage()
        sys.exit(1)
    
    myconf = Configuration(args.conf)
    DEfile = args.DEfile
    template_gff = args.template_gff
    if not DEfile or not template_gff:
        sys.stderr.write("Error: defile or template_gff is missing...\n")
        sys.exit(3)
    outpath = args.outpath
    godb = myconf.get("all").get("godb")
    kodb = myconf.get("all").get("kodb")
    func_annot(myconf,DEfile,template_gff,outpath=outpath,godb=godb,kodb=kodb)

if __name__ == '__main__':

    import sys

    main(sys.argv)
