#!/usr/bin/env python

import sys

'''when offered a map ID, this script can help to find the KO included in this map.
map2ko file is needed.
'''

map2kofile = open(sys.argv[1])
mapidfile = sys.stdin
warning = False
mapids = []
for eachline in mapidfile:
	tmp = eachline.strip().split()
	mapids.extend(tmp)

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
		mylist1.append(tmp[col1])
		mylist2.append(tmp[col2])
	return makedic_frm_list(mylist1,mylist2)
	
def get_koid(mapid,map2kodic):

	mydata = map2kodic.get("path:"+mapid)
	if not mydata:
		mydata = map2kodic.get(mapid)
	return mydata

map2kodic = makedic_frm_file(0,1,map2kofile,False)	
mydata = []
for eachmapid in mapids:
	try:
		mydata.extend(get_koid(eachmapid,map2kodic))
	except:
		if warning:
			print("[Warning] %s have no result"%eachmapid)
		else:
			pass
	
if not mydata:
	print("[error] can not find mapIDs:%s"%" ".join(mapids))
	sys.exit()
	
print("\n".join(mydata))
