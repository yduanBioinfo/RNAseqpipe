#!/usr/bin/env python3

'''install RNAseqpipe confs'''

import sys, os
import pkg_resources
from RNAseqpipe.install.init_base_conf import fill_prog_path
from RNAseqpipe.install.copy_conf import copy_static
from RNAseqpipe.install.copy_test import copy_test_files

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(\
        formatter_class = argparse.RawDescriptionHelpFormatter,\
        description = __doc__)
    parser.add_argument('-o','--outdir',nargs='?',help="output file",default='./')
    args = parser.parse_args(argv[1:])

    # Init base.conf
    conf_dir = pkg_resources.resource_filename('RNAseqpipe', 'confs/')
    base_template = pkg_resources.resource_filename('RNAseqpipe', 'confs/base.template')
    base_conf = pkg_resources.resource_filename('RNAseqpipe', 'confs/base.conf')
    with open(base_conf,'w') as out_conf:
        fill_prog_path(base_template, out_conf)

    # Copy confs and groups to dest directory
    copy_static(args.outdir)
    copy_test_files(args.outdir)

if __name__ == '__main__':

    main(sys.argv)

