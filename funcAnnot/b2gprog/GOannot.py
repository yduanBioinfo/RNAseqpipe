#!/usr/bin/env python

import sys, os
from RNAseqpipe.progsuit import Configuration, Prog_Rsp

def get_filename(path):
	
	return os.path.splitext(os.path.split(path)[1])[0]

def blast2go(conf,infile,outpath,gff=None,has_header=True,genome=None,silence=False):
	#infile: first column is ID[a str refers to file path]

	progname = "blast2go"
	order1 = {"--infile":infile}
	if not has_header:
		order1["-n"] = ""
	if gff:
		order1["-g"] = gff
	if genome:
		order1["-s"] = genome
	dbs = conf.get("pGOdb").values()
	order2 = {"-d":" ".join(dbs),"-o":outpath}
	
	blast2go = Prog_Rsp(conf,progname,order1,order2,silence) 
	blast2go.run()
	
	outpath = os.path.abspath(outpath)
	if not outpath.endswith("/"):
		outpath = outpath+"/"
	annotalldb = outpath+get_filename(infile)+"_all.annot"
	annotdb1 = outpath+get_filename(infile)+get_filename(dbs[0])
	return annotalldb, annotdb1
	
