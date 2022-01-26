# Test funcAnnot

import sys, os
import subprocess
import pytest
from ..test.function_for_test import check_exist

test_go_enrich_res = "./test/RNAseqpip_data/data/test_GO_enrich.txt"

def test_enrich_Bubble_base(tmpdir):
    outfile = tmpdir.join("test.pdf")
    order = "./funcAnnot/go/enrich_bubble.R {} -o {}".format(test_go_enrich_res, outfile).split()
    subprocess.call(order)
    check_exist(str(outfile),100)

def test_enrich_Bubble_stdin(tmpdir):
    outfile = tmpdir.join("test.pdf")
    order = "cat {} | ./funcAnnot/go/enrich_bubble.R - -o {}".format(test_go_enrich_res, outfile)
    subprocess.call(order,shell=True)
    check_exist(str(outfile),100)

def test_enrich_Bubble_resize(tmpdir):
    outfile = tmpdir.join("test.pdf")
    order = "./funcAnnot/go/enrich_bubble.R {} -o {} --width 20 --height 10 --units cm".format(test_go_enrich_res, outfile).split()
    subprocess.call(order)
    check_exist(str(outfile),100)
