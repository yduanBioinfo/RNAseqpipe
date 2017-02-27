#!/usr/bin/env python 

import sys, re

'''
_version=0.0.2
generate file used in R for kegg enrichment alanysis
count to pathway
./keggCountal.py growth_05 ./mydata/aaa pathway.ko outfile
./keggCountal.py [kegg annot file(mygenes)] [kegg annot file(all genes)] \
[path2KO file] [out file]
'''
mygenefile = open(sys.argv[1])
myuniversefile = open(sys.argv[2])
map2kofile = open(sys.argv[3])
try:
	outfile = open(sys.argv[4],'w')
except:
	outfile = sys.stdout

koPartten = re.compile(r'ko:(.*)')
#mapPartten = re.compile(r'path:(ko.*)')  a
#mapPartten = re.compile(r'path:(map.*)')  b
#mapPartten = re.compile(r'path:(.*)')#should write as a|b
mapPartten = re.compile(r'path:(ko.*|map.*)') 

mapData = {}
#mapData = {
#map00010:[ko111,ko222,ko333,...],
#map00020:[ko123,ko234,ko345,...],
#[],[],[],...}
countData = {}
#countData = {
#map:[count of all,count of degs]}

for each in map2kofile:#generate mapData
	tmp = each.strip().split()
	try:
		mypath = mapPartten.match(tmp[0]).group(1)
		myko = koPartten.match(tmp[1]).group(1)
	except:
		continue
	mapData.setdefault(mypath,[]).append(myko)

ko2mapData = {}
#ko2mapData = {
#ko123:[map00010,map233,map555,...],
#ko234:[map4556,map55666,map234,...].
#[],[],[]...}

for eachkey in mapData.keys():
	for eachv in mapData[eachkey]:
		ko2mapData.setdefault(eachv,[]).append(eachkey)

#get count of DEGs#
myko = ''
alofDEGs = 0
for each in mygenefile:
	try:
		myko = each.split()[1]
	except:
		continue
	if not ko2mapData.get(myko):
		#print('%s not found'%myko)
		continue
	for i in ko2mapData.get(myko):
		countData.setdefault(i,[0,0])[1] += 1
		alofDEGs += 1

#get count of all#
myko = ''
alofUniverse = 0
for each in myuniversefile:
	try:
		myko = koPartten.match(each.strip().split()[1]).group(1)
	except:
		try:
			myko = each.split()[1]
		except:
			continue

	if not ko2mapData.get(myko):
		#print
		continue
	for i in ko2mapData.get(myko):
		countData.setdefault(i,[0,0])[0] += 1
		alofUniverse += 1

def writefile(data):

	outfile.write("The map ID\tThe Counts of All\tNum of All\tThe Counts of DEGs\tNum of Al DEG\n")
	mykeys = data.keys()
	mykeys.sort()
	for i in mykeys:
		outfile.write("%s\t%d\t%d\t%d\t%d\n"%\
(i,data[i][0],alofUniverse,data[i][1],alofDEGs))

writefile(countData)
mygenefile.close()
myuniversefile.close()
map2kofile.close()
outfile.close()
