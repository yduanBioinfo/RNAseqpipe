#!/usr/bin/env python3

import sys
import re

blastout_start = re.compile(r"\s*\<BlastOutput_iterations\>\s*")
blastout_end = re.compile(r"(.*)\<\/BlastOutput_iterations\>\s*")
iteration_start = re.compile(r"\s*(\<Iteration\>)\s*")
iteration_end = re.compile(r"(\<\/Iteration\>).*")
def get_header(infile):
    header = []
    for line in infile:
        header.append(line)
        if blastout_start.match(line):
            return header

def iter_iteration(infile):

    def _first_line():
        iteration = []
        line = next(infile)
        if iteration_end.match(line):
            raise ValueError()
        #if not iteration_start.match(line):
        #    iteration.append(line)
        iteration.append(line)
        return iteration

    start = "<Iteration>\n"
    iteration = _first_line()
    for line in infile:
        ### End of blastout iteration
        _m = blastout_end.match(line)
        if _m:
            iteration.append(_m.group(1)+"\n")
            yield iteration
            return
        ###
        # End iteration
        _m = iteration_end.match(line)
        if _m:
            iteration.append(_m.group(1)+"\n")
            yield iteration
            iteration= _first_line()
            ### End of of blastout iteration
            if blastout_end.match(iteration[0]):
                return
            ###
            continue
        iteration.append(line)

def iter_body(infile,max_iter):
    data = []
    i = 0
    for block in iter_iteration(infile):
        ## aa
        data.append(block)
        i += 1
        if i == max_iter:
            yield data
            data = []
            i = 0
    yield data

def write_file(outfile,header,body):
    outfile.write("".join(header))
    for i in body:
        outfile.write("".join(i))
    return outfile

"""
Split blast out xml file into short files.
"""
infile = open(sys.argv[1])
outprefix = sys.argv[2]
## aa
max_iter = 50
#max_iter = 1000
#max_iter = 1000000000
header = get_header(infile)
# iterate body
outfiles = []
x=0
for body in iter_body(infile,max_iter):
    outfile = open(outprefix+"_"+str(x)+"_.xml",'w')
    outfiles.append(outfile)
    ## aa
    write_file(outfile,header,body)
    x += 1


## get tail
tail = ["</BlastOutput_iterations>\n"]
for line in infile:
    tail.append(line)

# output tails
for outfile in outfiles:
    outfile.write("".join(tail))
    outfile.close()
