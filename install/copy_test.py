#!/usr/bin/env python3

""" Copy test files to dest directory. """
import sys, os.path
import pkg_resources, shutil

def _copy_file(filename, dest):
    src = pkg_resources.resource_filename('RNAseqpipe', 'test/RNAseqpip_data/'+filename)
    shutil.copy(src, dest)

def copy_test_files(dest):
    out_test_dir_name = 'test'
    # Copy test data
    data_dir = pkg_resources.resource_filename('RNAseqpipe', 'test/RNAseqpip_data/data/')
    shutil.copytree(data_dir, os.path.join(dest, out_test_dir_name, 'test_bash', 'data'))
    
    # Copy scripts
    scripts = ['cmd', 'init_test_data.sh', 'run_test.sh', 'init_test_conf.sh']
    script_dest = os.path.join(dest, out_test_dir_name, 'test_bash')
    for i in scripts:
        _copy_file(i, script_dest)
