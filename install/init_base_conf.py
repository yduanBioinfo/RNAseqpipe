#!/usr/bin/env python3

from RNAseqpipe.progsuit import Configuration 
from shutil import which

def fill_prog_path(input_file, out_file):
    conf = Configuration(input_file)
    for key, val in conf.items():
        out_file.write("<%s>\n" % key)
        if key == 'prog_path':
            for subkey, subval in val.items():
                exc_path = which(subkey)
                if exc_path == None:
                    out_file.write("%s=%s\n" % (subkey, ""))
                else:
                    out_file.write("%s=%s\n" % (subkey, exc_path))
        else:
            for subkey, subval in val.items():
                out_file.write("%s=%s\n" % (subkey, subval))
        out_file.write('</%s>\n' % key)

if __name__ == '__main__':
    import sys
    fill_prog_path(sys.argv[1],sys.stdout)
