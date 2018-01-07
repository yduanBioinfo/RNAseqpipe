#!/usr/bin/env python

import sys

'''when offered a map ID, this script can help to find the KO included in this map.
map2ko file is needed.
'''
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

def map2ko(mapidfile,map2kofile,outfile,sep="\t",warning=False,keep_map_id=False):

    mapids = []
    for eachline in mapidfile:
        tmp = eachline.strip().split()
        mapids.extend(tmp)

    map2kodic = makedic_frm_file(0,1,map2kofile,False)    
    mydata = []
    def mode1():
    #do not keep mapid information
    #only output KO ids
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
    
        outfile.write("\n".join(mydata))
        outfile.write("\n")

    def mode2():
    #output mapid+sep+ko
        for eachmapid in mapids:
            kos = get_koid(eachmapid,map2kodic)
            if not kos:
                continue
            for eachko in kos:
                outfile.write(eachmapid+sep+eachko+"\n")
    if keep_map_id:
        mode2()
    else:
        mode1()

def main(argv):

    import argparse

    parser = argparse.ArgumentParser(description="[Usage] echo \"map00010\" | ./map2ko.py - -d ~/data/funcdb/kodb/pathway.ko")
    parser.add_argument('infile',nargs='?',help="file to be filtered, \"-\" for stdin ")
    parser.add_argument('-d','--map2ko',nargs='?',help="map2ko dbfile",type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('--keep',help="keep map id",default=False,action='store_true')
    #parser.add_argument('-s','--sep',nargs='?',default="\t")
    #parser.add_argument('-n1','--has_header1',action='store_false',help="if infile hasn't header,set this option(deprecate)")
    #parser.add_argument('-n2','--has_header2',action='store_false',help='same with n1 but for file2(deprecate)')
    #parser.add_argument('-m','--method',nargs='?',default='filter',choices=['filter','differ'])
    args = parser.parse_args(argv[1:])

#    sep1=args.sep1
#    sep2=args.sep2
#    keycol1=args.keycol1#key column which define filter limits
#    keycol2=args.keycol2
#    has_header1=args.has_header1
#    has_header2=args.has_header2
    
    if args.infile == '-':
        infile=sys.stdin
    else:
        infile=open(args.infile)

    map2ko(infile,args.map2ko,args.outfile,keep_map_id=args.keep)

if __name__ == '__main__':

    import sys
    main(sys.argv)
