import subprocess
import os
import errno

import pytest

def application_installed(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

@pytest.mark.dependencies
class TestClass:

    def test_matplotlib_installed(self):
        try:
            import matplotlib
        except ImportError:
            assert False, "matplotlib is not installed. You will not be able to use some commands."

    def test_bowtie_installed(self):
        assert application_installed('bowtie2'), "bowtie2 could not be found, and is needed for mapping raw sequencing data"

    def test_fastqc_installed(self):
        assert application_installed('fastqc'), "fastqc could not be found, and is needed for performing sample quality control"

    def test_fastq_screen_installed(self):
        assert application_installed('fastq_screen'), "fastq_screen could not be found, and is needed for performing sample quality control"

    def test_samtools_installed(self):
        assert application_installed('samtools'), "samtools could not be found, and is needed for removing PCR duplicates"

    def test_bedtools_installed(self):
        assert application_installed('bedtools'), "bedtools could not be found, and is needed for calling positive windows"

    def test_bedGraphToBigWig_installed(self):
        assert application_installed('bedGraphToBigWig'), "bedGraphToBigWig could not be found, and is needed for creating bigwig files"

    def test_bedToBigBed_installed(self):
        assert application_installed('bedToBigBed'), "bedToBigBed could not be found, and is needed for creating bigbed files"

