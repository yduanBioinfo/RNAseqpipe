# Test funcAnnot

import sys, os
import subprocess
import pytest
#from ..test.function_for_test import _check_exist

test_go_enrich_res = "./test/RNAseqpip_data/data/test_GO_enrich.txt"

def test_enrich_Bubble(tmpdir):
    outfile = tmpdir.join("test.pdf")
    order = "./funcAnnot/go/enrich_bubble.R {} {}".format(test_go_enrich_res, outfile).split()
    subprocess.call(order)
    assert os.path.exists(outfile) and os.path.getsize(outfile) > 0
