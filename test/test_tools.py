# Test funcAnnot

import sys, os
import subprocess
import pytest
import filecmp
from ..test.function_for_test import check_exist

test_input_gc_gtf = "./test/RNAseqpip_data/data/test.dlmrna.gc.final.gtf"
test_input_hs_gtf = "./test/RNAseqpip_data/data/test_human.gtf"

def test_get_gene_length_gc(tmpdir):
    ref_file = './test/RNAseqpip_data/refout/tools/test_get_gene_length_gc.txt'
    outfile = tmpdir.join("test.txt")
    order = "./get_gene_length.py {} -o {}".format(test_input_gc_gtf, outfile).split()
    subprocess.call(order)
    assert filecmp.cmp(outfile, ref_file)

def test_get_gene_length_hs(tmpdir):
    ref_file = './test/RNAseqpip_data/refout/tools/test_get_gene_length_hs.txt'
    outfile = tmpdir.join("test.txt")
    order = "./get_gene_length.py {} -o {}".format(test_input_hs_gtf, outfile).split()
    subprocess.call(order)
    assert filecmp.cmp(outfile, ref_file)

def test_get_gene_length_hs_tx(tmpdir):
    ref_file = './test/RNAseqpip_data/refout/tools/test_get_gene_length_hs_tx.txt'
    outfile = tmpdir.join("test.txt")
    order = "./get_gene_length.py {} -o {} --transcripts".format(test_input_hs_gtf, outfile).split()
    subprocess.call(order)
    assert filecmp.cmp(outfile, ref_file)

def test_get_gene_length_hs_featureName(tmpdir):
    ref_file = './test/RNAseqpip_data/refout/tools/test_get_gene_length_hs_fn.txt'
    outfile = tmpdir.join("test.txt")
    order = "./get_gene_length.py {} -t {} -o {}".format(test_input_hs_gtf, "db_xref", outfile).split()
    subprocess.call(order)
    assert filecmp.cmp(outfile, ref_file)
