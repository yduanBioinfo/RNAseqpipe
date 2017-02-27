#!/usr/bin/env python

import os

'''
_version=0.0.1
needed script:
subset_fasta.py
find_genes_in_genome.py
'''
def get_supp(list1,list2):#supplementary list

	outlist = list1[:]
	for each in list2:
		outlist.remove(each)
	return outlist
	
def get_filename(path):
	
	return os.path.splitext(os.path.split(path)[1])[0]
	
def subset_fasta(infasta,IDfile,outfasta):

	os.system("../subset_fasta.py %s %s %s"%(infasta,IDfile,outfasta))
	
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
		IDfile = open(outpath+get_filename(file)+".ID",'w')
		IDfile.write("\n".join(IDs))
		IDfile.close()
		
	return IDs
	
def get_sig_fasta(genome,gtf_file,IDfile,outfile):

	os.system("../find_genes_in_genome.py "+genome+" "+gtf_file+" "+IDfile+">"+outfile)
	
def blastx(genome,fasta,outfile):

	os.system("blastx -query "+fasta+" -db "+genome+" -out "+outfile+" -outfmt 5 -evalue 0.000001 -max_target_seqs 5 -num_threads 8")
	
def b2g(xmlfile,outfile):

	os.system(r"b2g -in "+xmlfile+" -annot -out "+outfile)
	
def fasta2annot(fastafile,db,allIDs,outpath=""):

	'''
	fasta -> xml -> annot -> IDs -> leftIDs -> unannoted fasta (-> xml)
	xml name:fastafile+db+.xml
	annot file name: fastafile+db+.annot
	leftIDs name: fastafile+db+.lftID
	unannoted fasta name: fastafile+db+_lft.fasta
	'''
	
	xmlfile = outpath+get_filename(fastafile)+get_filename(db)+".xml"
	blastx(db,fastafile,xmlfile)
	annotfile = xmlfile[:-4]#del ".xml" 
	b2g(xmlfile,annotfile)#blast2go out annotfile+".annot"
	annotedIDs = get_IDs(annotfile+".annot")
	leftIDs = get_supp(allIDs,annotedIDs)
	tmp = open(annotfile+".lftID",'w')#write left IDs
	tmp.write("\n".join(leftIDs))
	tmp.close()#write left IDs
	subset_fasta(fastafile,annotfile+".lftID",annotfile+".lft_fasta")
	return leftIDs
	
def main(argv):
	
	import argparse, sys, os
	
	parser = argparse.ArgumentParser(description='blast2go pip')
	parser.add_argument('--infile',help='differential expression annalys out file,\
	typically: NOISeq out file(the first column is gene_id and has header line)',\
	nargs='?')
	parser.add_argument('-n','--has_header',help='if infile hasn\'t header, set this option',\
	default=True,action='store_false')
	parser.add_argument('-g','--gff_file',nargs='?',help='gff file')
	parser.add_argument('-s','--genome',nargs='?',help='genome fasta file')
	parser.add_argument('-d','--dbs',nargs='+',help='protein dbs, .faa file')
	parser.add_argument('-o','--outpath',nargs='?',help='out file path',default='./b2g')
	args = parser.parse_args(argv[1:])	
	
	dbnames = map(get_filename,args.dbs)
	
	outpath = args.outpath
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

	owd = os.getcwd()#work directory when script initial
	#may need fix
	sdir = os.path.dirname(os.path.abspath(argv[0]))#script directory
	os.chdir(sdir)#can cause bug!!!
	
	myIDs = get_IDs(args.infile,args.has_header,writefile=True,outpath=outpath)#write ID file myIDs include all genes
	IDfile = outpath+get_filename(args.infile)+".ID" 
	
	sig_fasta = outpath+get_filename(args.infile)+".fasta"
	get_sig_fasta(args.genome,args.gff_file,IDfile,sig_fasta)#write fasta
	
	lftIDs = fasta2annot(sig_fasta,args.dbs[0],myIDs,outpath)
	fastaname = args.infile
	for i in range(1,len(args.dbs)):
		fastaname = outpath+get_filename(fastaname)+get_filename(args.dbs[i-1])+".lft_fasta"
		lftIDs = fasta2annot(fastaname,args.dbs[i],lftIDs,outpath)
	
	os.chdir(outpath)
	os.system(r"cat *.annot > %s_all.annot" % get_filename(args.infile))
	os.chdir(owd)#change work directory back
	
if __name__ == '__main__':

	import sys
	
	main(sys.argv)
	
