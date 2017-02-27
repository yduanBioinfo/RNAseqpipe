#!/usr/bin/env python

import sys

annotfile = open(sys.argv[1])
kofile = sys.stdin
koids = []
warning = False
for eachline in kofile:
	tmp = eachline.strip().split()
	koids.extend(tmp)
	
def makedic_frm_list(list1,list2):

	if len(list1) != len(list2):
		print("length error")
		sys.exit()
	mydic = {}
	for i in range(len(list1)):
		mydic.setdefault(list1[i],[]).append(list2[i])
	return mydic
	
def makedic_frm_file(col1,col2,file,has_header=True,sep="\t"):

	if has_header:
		file.readline()
	mylist1 = []
	mylist2 = []
	for eachline in file:
		tmp = eachline.strip().split(sep)
		try:#some line only have one element while not annoted
			mylist1.append(tmp[col1])
			mylist2.append(tmp[col2])
		except:
			pass
	return makedic_frm_list(mylist1,mylist2)
	
def get_geneid(geneid,ko2genedic):

	mydata = ko2genedic.get(geneid)
	if not mydata:
		mydata = ko2genedic.get(geneid.lstrip("ko:"))
	return mydata
	
ko2genedic = makedic_frm_file(1,0,annotfile,False)
mydata = []
for eachko in koids:
	try:
		mydata.extend(get_geneid(eachko,ko2genedic))
	except:
		if warning:
			print("[Warning] %s have no result"%eachko)
		else:
			pass
	
if not mydata:
	print("[error] can not find ko(s):%s"%" ".join(koids))
	sys.exit()

print("\n".join(mydata))