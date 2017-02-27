#!/usr/bin/env python

'''
_version = 0.2.1
'''

import sys, os, copy
from collections import OrderedDict as Ordic
from progsuit import Configuration, Prog_Rsp

def cufflink(conf,bamf,outpath,silence=False,thre=0):
	#conf if configuration obj.
	#bamf: bam file [path]
	#thre : -p value
	progname = "cufflinks"
	order1 = {bamf:"","-o":outpath}
	order2 = {}
	cufflink = Prog_Rsp(conf,progname,order1,order2,silence)
	return cufflink.run()

def cufflinks(conf,bamfs,outpath=None,silence=False,maxp=20):
	#conf is a Configuration obj
	#bamfs: bam/sam files [list]
	#if outpath there isn't a specific outpath,
	#out transcripts.gtf should write in the same dir of bam file.
	import time
	assert isinstance(bamfs,list)
	
	def child():
		cufflink(myconf,bamf,path,silence)
		os._exit(0)
		
	myconf = conf
	outfiles = []
	if outpath and len(bamfs) != len(outpath):
		error_handle(1,"cufflinks")
	
	try:
		thre = myconf["cufflinks"]["-p"]#thread for cufflink: -p
	except:
		thre = 4
	max_ = maxp // int(thre)#max_: max subprocess should be created
	pids = []#save pid of subprocess
	for i in range(len(bamfs)):
		bamf = bamfs[i]
		try:
			path = outpath[i]
		except:
			path = os.path.dirname(bamf)
		
		pid = os.fork()
		if pid <= 0:
			break
		outfiles.append(path+"/transcripts.gtf")#parent process
		
		pids.append(pid)
		if len(pids) == max_:
			myid, status = os.wait()
			pids.remove(myid)
			
	if pid == 0:
		child()
		
	try:#waiting for all threads done
		while 1:
			os.wait()
	except OSError:
		pass
		
	return outfiles

def cuffmerge(conf,assembly,outpath,silence=False):
	#assembly: assembly.txt
	
	order1 = Ordic([("-o",outpath)])
	order2 = Ordic()
	if conf["all"]["gff"]:
		order2["-g"] = conf["all"]["gff"]
	if conf["all"]["genome"]:
		order2["-s"] = conf["all"]["genome"]
	order2[assembly] = ""
	cuffmerge = Prog_Rsp(conf,"cuffmerge",order1,order2,silence)
	cuffmerge.run()
	del cuffmerge
	return outpath+"/merged.gtf"
	
def cuffquant(conf,bamf,gtffile,outpath,silence=False):

	order1 = {"-o":outpath,"locp":[gtffile,bamf]}
	order2 = {"-p":8}
	Pcuffquant = Prog_Rsp(conf,"cuffquant",order1,order2,silence)
	Pcuffquant.run()

def cuffnorm(conf,bamfs,gtffile,outpath,silence=False):
	#bamfs can be .bam .sam .cxb 
	
	assert isinstance(bamfs,list)
	outpath = outpath+"/norm_out"
	mylocp = [gtffile]
	mylocp.extend(bamfs)
	order1 = {"-o":outpath,"locp":mylocp}
	order2 = {"-p":8}
	Pcuffnorm = Prog_Rsp(conf,"cuffnorm",order1,order2,silence)
	Pcuffnorm.run()
	return outpath
	
def cuffquants(conf,bamfs,gtffile,outpath=None,silence=False):

	assert isinstance(bamfs,list)
	
	if outpath and len(bamfs) != len(outpath):
		error_handle(1,"cuffquants")
	
	outfiles = []
	for i in range(len(bamfs)):
		bamf = bamfs[i]
		try:
			path = outpath[i]
		except:
			path = os.path.dirname(bamf)
		cuffquant(conf,bamf,gtffile,path,silence)
		outfiles.append(path+"/abundances.cxb")
	return outfiles

def cuffdiff(conf,gff,bams1,bams2,outpath=None,silence=False):

	if not outpath:
		outpath = "./"
	order1 = [(gff,""),(",".join(bams1),""),(",".join(bams2),""),("-o",outpath)]
	order2 = {"-p":8}
	progname = "cuffdiff"
	prog = Prog_Rsp(conf,progname,order1,order2,silence)
	prog.run()
	
def error_handle(code,name="unknown"):

	if code == 1:
		sys.stderr.write("Error: outpaths should be as many as bam/sam files are")
	if code:
		sys.stderr.write("Error in %s step.\n" % name)
		sys.exit(20)

if __name__ == '__main__':			
	
	myconf = Configuration(open("../conf2.txt"))
	bam = "/home/yduan/data/growth/ruibo/pip/tmp/2014-11-BB2/sort.bam"
	cufflink(myconf,bam,"aaaa")