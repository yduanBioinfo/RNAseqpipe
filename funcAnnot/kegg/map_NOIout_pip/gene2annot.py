#!/usr/bin/env python

import sys

myexpfile = open(sys.argv[1])
geneidfile = sys.stdin
geneids = []
for eachline in geneidfile:
	tmp = eachline.strip().split()
	geneids.extend(tmp)

def makedic_frm_list(list1,list2):

	if len(list1) != len(list2):
		print("length error")
		sys.exit()
	mydic = {}
	for i in range(len(list1)):
		mydic.setdefault(list1[i],[]).append(list2[i])
	return mydic
	
def makedic_frm_exp(file,sep="\t"):#return a dic of expression file(NOIout)

	firstline = file.readline()
	mylist1 = []
	mylist2 = []
	for eachline in file:
		tmp = eachline.strip().split(sep)
		mylist1.append(tmp[0])
		mylist2.append(sep.join(tmp[1:]))
	return firstline, makedic_frm_list(mylist1,mylist2)
	
def write_exp(expftln,expdic,geneids):

	print(expftln.strip())
	for eachid in geneids:
		tmp_exp = expdic.get(eachid)
		if tmp_exp:
			print(eachid+"\t"+tmp_exp[0])
		
expftln, expdic = makedic_frm_exp(myexpfile)
write_exp(expftln,expdic,geneids)