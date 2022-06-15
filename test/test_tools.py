# Test funcAnnot

import sys, os
import subprocess
import pytest
from ..test.function_for_test import check_exist

test_input_gc_gtf = "./test/RNAseqpip_data/data/test.dlmrna.gc.final.gtf"
test_input_hs_gtf = "./test/RNAseqpip_data/data/test_human.gtf"

def test_get_gene_length_gc(tmpdir):
    outfile = tmpdir.join("test.txt")
    order = "./get_gene_length.py {} -o {}".format(test_input_gc_gtf, outfile).split()
    subprocess.call(order)
    check_exist(str(outfile),100)

def test_get_gene_length_hs(tmpdir):
    outfile = tmpdir.join("test.txt")
    order = "./get_gene_length.py {} -o {}".format(test_input_hs_gtf, outfile).split()
    subprocess.call(order)
    check_exist(str(outfile),100)

def test_get_gene_length_hs_tx(tmpdir):
    outfile = tmpdir.join("test.txt")
    order = "./get_gene_length.py {} -o {} --transcripts".format(test_input_hs_gtf, outfile).split()
    subprocess.call(order)
    check_exist(str(outfile),100)
