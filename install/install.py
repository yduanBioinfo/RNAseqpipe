#!/usr/bin/env python3

'''install RNAseqpipe confs'''

import sys, os
import pkg_resources
from progsuit import log
from RNAseqpipe.install.init_base_conf import fill_prog_path
from RNAseqpipe.install.copy_conf import copy_static
from RNAseqpipe.install.copy_test import copy_test_files


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__)
    parser.add_argument('-o', '--outdir', nargs='?',
                        help="output file", default='./')
    parser.add_argument('-u', '--user-only',
                        help="Only copy and init for user",
                        action='store_true')
    args = parser.parse_args(argv[1:])

    if not args.user_only:
        print("Run step A")
        # Init base.conf
        base_template = pkg_resources.resource_filename('RNAseqpipe',
                                                        'confs/base.template')
        base_conf = pkg_resources.resource_filename('RNAseqpipe',
                                                    'confs/base.conf')
        try:
            with open(base_conf, 'w') as out_conf:
                fill_prog_path(base_template, out_conf)
        except PermissionError:
            log.warning("System level Base conf file was not changed"
                        "due to permission limits.")

    # Copy confs and groups to dest directory
    copy_static(args.outdir)
    copy_test_files(args.outdir)

    print("Run step B")
    # Init local base files.
    base_template = os.path.join(args.outdir, 'confs/base.template')
    base_conf = os.path.join(args.outdir, 'confs/base.conf')
    with open(base_conf, 'w') as out_conf:
        fill_prog_path(base_template, out_conf)


if __name__ == '__main__':

    main(sys.argv)
