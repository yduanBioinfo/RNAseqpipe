#!/usr/bin/env python3

"""
  Merge flag stat results of samtools

  [Usage]: ./merge_flagestat.py flagstat1.txt flagstat2.txt ... > merged_flagstat.txt

  input: Series of flagstat.txt file
    98692501 + 0 in total (QC-passed reads + QC-failed reads)
    2482577 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    90171005 + 0 mapped (91.37% : N/A)
    96209924 + 0 paired in sequencing
    48104962 + 0 read1
    48104962 + 0 read2
    84003548 + 0 properly paired (87.31% : N/A)
    84863732 + 0 with itself and mate mapped
    2824696 + 0 singletons (2.94% : N/A)
    589024 + 0 with mate mapped to a different chr
    525005 + 0 with mate mapped to a different chr (mapQ>=5)
  output:
    # QC-passed reads
    Metics	File1	File2	File3
    Total	98692501	98377103	87672384
    Secondary	...	...	...
    ...	...	...	...
"""

import sys, re
from collections import OrderedDict as ordic

class Flagstat(ordic):

    """
        Read flagstat file, and saved into dict-like data.
        Data:{key:(QCpassed,QCunpassed)}
    """

    sep = r"(?:\b|\s+)"
    line_with_bracket = re.compile(sep.join(["^(\d+)","\+","(\d+)\s+([^(]*)","\((\d+\.\d+%)",":","([^:]+)\)$"]))
    line_pattern = re.compile(sep.join([r"^(\d+)","\+","(\d+)\s+([^(]+)",".*"]))

    def __init__(self,infile):
        for line in infile:
            rbp = Flagstat.line_with_bracket.match(line)
            rp = Flagstat.line_pattern.match(line)
            if rbp:
                self._save_rbp(rbp)
            elif rp:
                self._save_rp(rp)
            else:
                raise KeyError("This line is not support:%s\n" % line)

    # line with bracket.
    # 90171005 + 0 mapped (91.37% : N/A)
    def _save_rbp(self,rp):
        self._save_rp(rp)
        # Save rates.
        self[rp.group(3)+"_rate"] = (rp.group(4),rp.group(5))

    def _save_rp(self,rp):
        self[rp.group(3)] = (rp.group(1),rp.group(2))

# Guess sample id from path.
def guess_id(name,idx=-2):
    name_array = name.strip().split("/")
    return name_array[idx]

# Only output passed qc data
# stats: list of Flagstat objects.
# header: list of filename
def write_QC(stats,outfile,header=[],sep="\t"):
    if header:
        outfile.write("Metric"+sep+sep.join(header)+"\n")
    for key in stats[0].keys():
        # x.get(key)[0]: passed qc; [1] unpassed.
        outfile.write(key+sep+sep.join(map(lambda x:x.get(key)[0],stats))+"\n")
infiles = sys.argv[1:]
ffs = list(map(lambda x:Flagstat(open(x)),infiles))
names = list(map(guess_id,infiles))

write_QC(ffs,sys.stdout,header=names)
