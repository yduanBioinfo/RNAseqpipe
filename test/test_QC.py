# Test funcAnnot

import sys, os
import subprocess
import pytest
from ..test.function_for_test import check_exist, check_exists

"""../run_RNAseqpipe.py QC ../test/RNAseqpip_data/data/fqs/test1_1.fq ../test/RNAseqpip_data/data/fqs/test1_2.fq ../test/RNAseqpip_data/data/fqs/test2_1.fq ../test/RNAseqpip_data/data/fqs/test2_2.fq ../test/RNAseqpip_data/data/fqs/test3_1.fq ../test/RNAseqpip_data/data/fqs/test3_2.fq -o xxxx -t 10"""

def _gen_fastqcfile(infiles, num):
    outdata = []
    for i in infiles:
        # remove '.fq.gz' or '.fq'
        outdata.append(i[:-num]+'_fastqc.html')
        outdata.append(i[:-num]+'_fastqc.zip')
    return outdata

def gen_fastqcfile_fqgz(infiles):
    return _gen_fastqcfile(infiles, 6)

def gen_fastqcfile_fq(infiles):
    return _gen_fastqcfile(infiles, 3)

def test_QC(tmpdir):
    order = "./run_RNAseqpipe.py QC ./test/RNAseqpip_data/data/fqs/test1_1.fq ./test/RNAseqpip_data/data/fqs/test1_2.fq ./test/RNAseqpip_data/data/fqs/test2_1.fq ./test/RNAseqpip_data/data/fqs/test2_2.fq ./test/RNAseqpip_data/data/fqs/test3_1.fq ./test/RNAseqpip_data/data/fqs/test3_2.fq -o {} -t {}".format(tmpdir, 10).split()
    files = ['test1_1.fq','test1_2.fq','test2_1.fq','test2_2.fq','test3_1.fq','test3_2.fq']
    clean_files = ['test1_1.clean_read1.fq.gz', 'test2_1.clean_read1.fq.gz',  'test3_1.clean_read1.fq.gz', 'test1_2.clean_read2.fq.gz', 'test2_2.clean_read2.fq.gz', 'test3_2.clean_read2.fq.gz']
    subprocess.call(order)
    # Check clean files
    check_exists(map(lambda x:os.path.join(tmpdir,'clean',x), clean_files), 200)
    # Check reports
    check_exists(map(lambda x:os.path.join(tmpdir,'report_before_qc/fastqc',x), 
        gen_fastqcfile_fq(files)), 200)
    check_exists(map(lambda x:os.path.join(tmpdir,'report_after_qc/fastqc',x), 
        gen_fastqcfile_fqgz(clean_files)), 200)
    check_exist(os.path.join(tmpdir, 'report_before_qc/multiQC/multiqc_report.html'), 200)
    check_exist(os.path.join(tmpdir, 'report_after_qc/multiQC/multiqc_report.html'), 200)
