#!/usr/bin/env python3 

import sys, re
import argparse

'''
_version=0.0.2
generate file used in R for kegg enrichment alanysis
count to pathway
./keggCountal.py growth_05 ./mydata/aaa pathway.ko outfile
./keggCountal.py [kegg annot file(mygenes)] [kegg annot file(all genes)] \
[path2KO file] [out file]
'''
"""
ko	ko:KXXXX or KXXXX
map	path:koxxxxx or path:mapxxxxx or mapxxxxx

pathway.ko file:
... ...
... ...
bg file:
... ...
... ...
"""
koPartten = re.compile(r'(ko:)?([^:]*)')
mapPartten = re.compile(r'(path:)?(ko.*|map.*)') 
SEP = "\t"

def writefile(data,outfile,alofUniverse,alofDEGs):

    outfile.write("The map ID\tThe Counts of All\tNum of All\tThe Counts of DEGs\tNum of Al DEG\n")
    mykeys = list(data.keys())
    mykeys.sort()
    for i in mykeys:
        outfile.write("%s\t%d\t%d\t%d\t%d\n"% (i,data[i][0],alofUniverse,data[i][1],alofDEGs))

def iter_map2ko(map2kofile):
    for each in map2kofile:
        tmp = each.strip().split(SEP)
        try:
            mypath = mapPartten.match(tmp[0]).group(2)
            myko = koPartten.match(tmp[1]).group(2)
        except:
            continue
        yield mypath, myko

def iter_la(infile):
    for line in infile:
        yield line.strip().split(SEP)

def get_count_of_DEGs(mygenefile,ko2mapData,countData,gene2ko=None):
    """
        When gene2ko is provide, get ko from gene/tx ID.
        Otherwise get ko from the second field of eachline.
    """
    alofDEGs = 0
    for la in mygenefile:
        try:
            if gene2ko:
                myko = gene2ko[la[0]]
            else:
                myko = koPartten.match(la[1]).group(2)
        except:
            continue
        if not ko2mapData.get(myko):
            continue
        for i in ko2mapData.get(myko):
            countData.setdefault(i,[0,0])[1] += 1
            alofDEGs += 1
    return alofDEGs

def get_count_of_universe(myuniversefile,ko2mapData,countData):
    alofUniverse = 0
    gene2ko = {}
    for la in myuniversefile:
        try:
            myko = koPartten.match(la[1]).group(2)
            gene2ko[la[0]] = myko
        except:
            continue
    
        if not ko2mapData.get(myko):
            continue
        for i in ko2mapData.get(myko):
            countData.setdefault(i,[0,0])[0] += 1
            alofUniverse += 1
    return alofUniverse,gene2ko

def count(mygenefile,myuniversefile,map2kofile,outfile,mode):
    """
        mode g: get DEGs' ko from gene/tx ID.
        mode a: get DEGs' ko from ko annotation.
    """
    
    mapData = {}
    #mapData = {
    #map00010:[ko111,ko222,ko333,...],
    #map00020:[ko123,ko234,ko345,...],
    #[],[],[],...}
    countData = {}
    #countData = {
    #map:[count of all,count of degs]}
    ko2mapData = {}
    #ko2mapData = {
    #ko123:[map00010,map233,map555,...],
    #ko234:[map4556,map55666,map234,...].
    #[],[],[]...}
    
    # generate mapData
    #for mypath, myko in iter_map2ko(map2kofile):
    for mypath, myko in map2kofile:
        mapData.setdefault(mypath,[]).append(myko)
    
    # generate ko2mapData
    for eachkey in mapData.keys():
        for eachv in mapData[eachkey]:
            ko2mapData.setdefault(eachv,[]).append(eachkey)
    
    alofUniverse, gene2ko = get_count_of_universe(myuniversefile,ko2mapData,countData)
    if mode == 'g':
        _gene2ko = gene2ko
    else:
        _gene2ko = None
    alofDEGs = get_count_of_DEGs(mygenefile,ko2mapData,countData,_gene2ko)
    writefile(countData,outfile,alofUniverse,alofDEGs)

def main(argv):

    parser = argparse.ArgumentParser(description="Get count table of KEGG which is input file for KEGG enrichment.")
    parser.add_argument('infile',nargs='?',type=argparse.FileType('r'),help="KEGG annotation of DEGs,(or simply genelist of DEGs)")
    parser.add_argument('bg',nargs='?',type=argparse.FileType('r'),help="Back ground of KEGG annotation")
    parser.add_argument('pthko',nargs='?',type=argparse.FileType('r'),help="File records pathway to KO relationship.")
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-m','--mode',nargs='?',help="input file mode, annot(a) or genelist(g)",default='a',choices=['a','g'])
    args = parser.parse_args(argv[1:])
    # Use iter_la and iter_map2ko makes count avilable to list-like object.
    count(iter_la(args.infile),iter_la(args.bg),iter_map2ko(args.pthko),args.outfile,args.mode)

if __name__ == '__main__':
    main(sys.argv)


