#!/usr/bin/env python

import sys, os, time
from RNAseqpipe.progsuit import Configuration, Prog_Rsp

def get_filename(path):
	
	return os.path.splitext(os.path.split(path)[1])[0]
	
def get_IDs(file,header=False,col=0,writefile=False,outpath=""):

	'''file : NOISeq file [defualt]
	   open NOISeq file and write the first col(ID)
	'''
	
	if isinstance(file,str):
		NOIfile = open(file)
	else:
		NOIfile = file
	IDs = []
	if header:
		NOIfile.readline()
	for each in NOIfile:
		tmp = each.split("\t")
		if tmp and tmp[col].strip() not in IDs:
			IDs.append(tmp[col].strip())
	if writefile:
		outpath = outpath if outpath.endswith("/") else outpath+"/"
		os.mkdir(outpath)#16/8/16
		IDfile = open(outpath+get_filename(file)+".ID",'w')
		IDfile.write("\n".join(IDs))
		IDfile.close()
		
	return IDs

def get_sig_fasta(genome,gtf_file,IDfile,outfile):
	#the abspath should be fixed
	os.system("/home/yduan/script/python/RNAseqpip/funcAnnot/find_genes_in_genome.py "+genome+" "+gtf_file+" "+IDfile+">"+outfile)
	
def koannot(conf,infile,outpath,gff=None,has_header=True,genome=None,silence=False):

	myIDs = get_IDs(infile,has_header,writefile=True,outpath=outpath)#write ID file myIDs include all genes
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

	IDfile = outpath+get_filename(infile)+".ID" 
	#print("where id file: %s"%IDfile)
	#print("the outpath:%s"%outpath)
	sig_fasta = outpath+get_filename(infile)+".fasta"
	
	if not gff:
		gff=conf.get("all").get("gff")
	if not genome:
		genome=conf.get("all").get("genome")
	get_sig_fasta(genome,gff,IDfile,sig_fasta)#write fasta
	
	progname = "kaas_online"
	order1 = {sig_fasta:""}
	if not has_header:
		order1["-n"] = ""
	outfile = outpath+get_filename(infile)+".koannot"
	outfile = outfile+time.strftime("_%m_%d",time.localtime()) \
	if os.path.isfile(outfile) else outfile
	order2 = {"-o":outfile}
	
	prog = Prog_Rsp(conf,progname,order1,order2,silence) 
	prog.run()
	return outfile
