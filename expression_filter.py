#!/usr/bin/env python

def filt_exp(exps,cutoff):
    #set each in line to 0 when mean < cutoff

    assert isinstance(exps,list)
    lenth = len(exps)
    exps = map(float,exps)
    mysum = sum(exps)
    if mysum/lenth < cutoff:
        #filt out
        exps = [0 for i in range(lenth)]
        return True, exps
    return False, exps

def filter(infile,outfile,set_0=False,has_header=True,sep="\t",thr=1):
    #thr: threshold
    #--set_0: set all to 0 instead of del the lines which are fail to pass the cutoff

    if has_header:
        header = infile.readline()
    outfile.write(header)

    num_filt = 0
    for eachline in infile:
        tmp = eachline.strip().split(sep)
        exps = tmp[1:]
        is_filted, exps = filt_exp(exps,thr)
        if is_filted and not set_0:#this line should be filted
            num_filt += 1
        else:
            x = tmp[0:1]
            x.extend(map(str,exps))
            outfile.write(sep.join(x)+"\n")
            if is_filted:#and set_0: this line be set to 0
                num_filt += 1
    sys.stderr.write("%d lines been filted\n"%num_filt)

def main(argv):

    import argparse, sys

    parser = argparse.ArgumentParser(description='this script is a filter of expression file')
    parser.add_argument('infile',nargs='?',default=sys.stdin,type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-c','--cutoff',nargs='?',default=1,type=float)
    parser.add_argument('-s','--set_0',default=False,action='store_true',help='set this argument when you want to set all to 0 instead of del the lines which are fail to pass the cutoff')
    args = parser.parse_args(argv)

    filter(args.infile,args.outfile,set_0=args.set_0,thr=args.cutoff)

if __name__ == '__main__':

    import sys
    
    main(sys.argv[1:])
