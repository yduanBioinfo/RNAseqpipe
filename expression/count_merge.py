#!/usr/bin/env python

import os

def gen_names(files,mod=2):
    #file:[fileobj1,fileobj2,...]
    #mod: 
    #files:[/path/to/subp1/name1.count, /path/to/subp1/name1.count,...]
    #mod1: colname is name1,name2,....
    #mod2: colname is subp1,subp2,....
    
    assert isinstance(files,list)
    outlst = []
    for each in files:
        name = each.name
        if mod == 1:
            name = os.path.basename(name)
            name = name.split(".")[0]
        elif mod == 2:
            name = name.split("/")[-2]
        else:
            sys.stderr.write("mod value error!!!\n")
            raise ValueError
        outlst.append(name)
    return outlst
    
# lfile: length file
def writef(files,outfile,no_header=False,sep="\t",com="_",lfile=None):
    #Start charactor of comment. Do not write out.
    
    if no_header:
        #Don't write header out.
        #Read, but not write out.
        mytmp = list(map(lambda x:x.readline().split(),files))

    # merge lfile to count file.
    # lfile default without header.
    if lfile:
        files.append(lfile)

    mytmp = list(map(lambda x:x.readline().split(),files))
    mylength = len(mytmp)
    try:
        while mytmp:
            tmp = mytmp[0][0]
            mytmp = [i[1] for i in mytmp]
            mytmp.insert(0,tmp)
            if not tmp.startswith(com):#name not start with com
                outfile.write(sep.join(mytmp)+"\n")
            mytmp = list(map(lambda x:x.readline().split(),files))
    except:
        pass
        
def main(argv):

    import argparse, sys
    
    parser = argparse.ArgumentParser(description='merge htseq-count out file')
    parser.add_argument('countf',help='count files',nargs='+',type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',help='outfile default: stdout',\
    default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-l','--length_file',nargs='?',help='Gene length file',\
    type=argparse.FileType('r'))
    parser.add_argument('-e','--none_header',help='not to write header',default=False,\
    action='store_true')
    parser.add_argument('-f','--no_ori_header',help='Ori count file contains header, but don\'t write it out.',\
    default=False,action='store_true')
    parser.add_argument('-n','--names',help='sample names to write on header',nargs='+')
    args = parser.parse_args(argv)
    
    if args.names:
        assert isinstance(args.names,list)
        if len(args.names) > 1 and len(args.names) != len(args.countf):
            sys.stderr.write("Error: --names should be as many as --countf\n")
            raise TypeError
        header = args.names
    else:
        header = gen_names(args.countf)
    header.insert(0,"gene")
    sep = "\t"
    if args.length_file:
        header.append("length")
    if not args.none_header:#write header
        args.outfile.write(sep.join(header)+"\n")
    writef(args.countf,args.outfile,args.no_ori_header,lfile=args.length_file)
    
if __name__ == '__main__':

    import sys
    
    main(sys.argv[1:])
