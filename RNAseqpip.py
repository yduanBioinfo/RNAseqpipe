#!/usr/bin/env python

'''
_version = 0.2.3
'''
from __future__ import print_function
import sys, re, os, copy
from threading import Thread

import hisat, cufflinks, NOISeq, htseq
import funcAnnot.b2gprog.GOannot as GOannot
import funcAnnot.kegg.KOannot as KOannot
from funcAnnot.kegg.annot2go_stat import annot2go_stat
from funcAnnot.annot_from_db import annot
from funcAnnot.kegg.koenrich.keggenrich import enrich as koenrich
from progsuit import Configuration, Group_data
from get_gene_length import len_for_Rsp
from get_geneids import get_geneids

def getAbsPath(inpath,default='./tmp'):

	if not inpath:
		inpath = default
	if not os.path.isdir(inpath):
		os.makedirs(inpath)
	return os.path.abspath(inpath)

def matchpath(names1,names2,mypath):#hidden bugs
	#retrieve path for names1 or names2 in mypath(list)
	#mostly required by DE step.

	assert isinstance(names1,list) and isinstance(names2,list) and isinstance(mypath,list)
	
	def getmatch(names,paths):
		outpath = []
		for name in names:
			for path in paths:
				if re.search(r"/"+name+r"/",path):outpath.append(path)
		return outpath
		
	paths1 = getmatch(names1,mypath)
	paths2 = getmatch(names2,mypath)
	return paths1,paths2
	
def seq2exp(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff=True):
	#from fastq to assembly gff
	#run_cdiff: if cuffdiff step should be performed

	ali_ress,sort_ress = hisat.pip_hisats(myconf,myfq1,myfq2,fqnames,ali_path,ali_name)
	#alignment results,samtools sorted results

	cufflinks_ress = cufflinks.cufflinks(myconf,sort_ress)#cufflinks results

	#write path/assembly.txt
		
	assembly = ali_path+"/assembly.txt"
	cu_ress_f = open(assembly,'w')
	cu_ress_f.write("\n".join(cufflinks_ress))
	cu_ress_f.close()
		
	merged = cufflinks.cuffmerge(myconf,assembly,ali_path+"/merged_asm")
	quants = cufflinks.cuffquants(myconf,sort_ress,merged)
	norm = cufflinks.cuffnorm(myconf,quants,merged,ali_path)
	gene_fpkm = norm+"/genes.fpkm_table"#cufflinks fpkm
	samples_table = norm+"/samples.table"
	#run cuffdiff
	treat_group = mygroup_data.get_g_cdgp()
	#{group1:[dataname1,dataname2,...],group2:[],...}
	try:
		gpnames1,gpnames2=treat_group.values()
	except:
		print(treat_group.values(),file=sys.stderr)
		run_cdiff = False
		sys.stderr.write("Warning:only two group data can run cuffdiff!\n")
		sys.stderr.write("cuffdiff step skipped\n")
	
	if run_cdiff:
		quants1,quants2 = matchpath(gpnames1,gpnames2,quants)
		cu_diff = ali_path+"/diffout"
		cufflinks.cuffdiff(myconf,merged,quants1,quants2,cu_diff)
	
	#cuffdiff end
	htcount = htseq.htseqpip(myconf,sort_ress,ali_path,merged,plotqa=False)
	#return expressionf	
	return gene_fpkm, htcount, quants, merged
	
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
	length = len_for_Rsp(mergedfile)
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
	
if __name__ == '__main__':

	import argparse
	
	parser = argparse.ArgumentParser(description='RNA-seq analyse pip')
	parser.add_argument('program',help='all for whole pip/ali for alignment/ass for assembly/\
	cptDE for compute different expression',choices=['all','cptDE','seq2exp','func'])
	parser.add_argument('-c','--conf',help='configuration file',nargs='?',\
	type=argparse.FileType('r'),default='configuration.txt')
	parser.add_argument('-1','--fq1',help='fq_1',nargs='*')
	parser.add_argument('-2','--fq2',help='fq_2',nargs='*')
	parser.add_argument('-g','--group_data',help='group_data file. conflict with -1 -2'\
	,nargs='?',type=argparse.FileType('r'))
	parser.add_argument('-s','--samples_table',help='cuffnorm samples.table file'\
	,nargs='?')
	parser.add_argument('-e','--expressionf',help='when chosing cptDE, this option is required.'\
	,nargs='?')
	parser.add_argument('-d','--DEfile',help='DEfile(NOISeq out). when chosing func,\
	this option is required.',nargs='?')
	parser.add_argument('-t','--template_gff',help='when chosing cptDE or func,\
	this option is used.',nargs='?')
	parser.add_argument('-o','--outpath',help='outpath',nargs='?')
	args = parser.parse_args(sys.argv[1:])

#	samples_table = None
#	if args.samples_table:
#		samples_table = args.samples_table
	myconf = Configuration(args.conf)
	ali_path = getAbsPath(args.outpath)#home path for alignment results
	ali_name = 'mapped.sam'#alignment result name
	
	if args.group_data and args.fq1:
		sys.stderr.write("You can't specify both -g and -1 the same time")
		raise StandardError
		
	if args.group_data:
		mygroup_data = Group_data(args.group_data)
		fqnames, myfq1, myfq2 = mygroup_data.get2pip()
	else:
		myfq1 = args.fq1
		myfq2 = args.fq2
		fqnames = ["" for i in range(len(myfq1))]
		
	if args.fq2 and len(args.fq1) != len(args.fq2):
		sys.stderr.write("-2 should be as long as -1\n")
		sys.exit(2)

	if args.program == 'all':
		
		gene_fpkm, htcount, quants, merged = \
		seq2exp(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data)
		nullexp,upexp,downexp = \
		cptDE(myconf,ali_path,mygroup_data,htcount,None,merged)
		func_annot(myconf,merged)
		
	if args.program == 'cptDE':
		
		nullexp,upexp,downexp = cptDE(myconf,ali_path,mygroup_data,args.expressionf,args.samples_table,args.template_gff)
		
	if args.program == 'seq2exp':
	
		seq2exp(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff=False)
		
	if args.program == 'func':
	
		DEfile = args.DEfile
		template_gff = args.template_gff
		if not DEfile or not template_gff:
			sys.stderr.write("Error: defile or template_gff is missing...\n")
			sys.exit(3)
		outpath = args.outpath
		godb = myconf.get("all").get("godb")
		kodb = myconf.get("all").get("kodb")
		func_annot(myconf,DEfile,template_gff,outpath=outpath,godb=godb,kodb=kodb)
	