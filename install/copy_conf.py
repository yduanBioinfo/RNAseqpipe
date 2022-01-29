#!/usr/bin/env python3

""" Copy static files to dest directory. """
import sys, os.path
import pkg_resources, shutil

def copy_static(dest):
    conf_dir = pkg_resources.resource_filename('RNAseqpipe', 'confs/')
    group_dir = pkg_resources.resource_filename('RNAseqpipe', 'group_data/')
    shutil.copytree(conf_dir, os.path.join(dest,'confs'))
    shutil.copytree(group_dir, os.path.join(dest,'group_data'))
