import subprocess
import os

def application_installed(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def test_pybedtools_installed():
    try:
        import pybedtools
    except ImportError:
        assert False, "pybedtools is not installed. You will not be able to use the 'plot_np' command."

def test_matplotlib_installed():
    try:
        import matplotlib
    except ImportError:
        assert False, "matplotlib is not installed. You will not be able to use the 'plot_np' command."

def test_metaseq_installed():
    try:
        import metaseq
    except ImportError:
        assert False, "metaseq is not installed. You will not be able to use the 'plot_np' command."

def test_bowtie_installed():
    assert application_installed('bowtie2')

def test_fastqc_installed():
    assert application_installed('fastqc'), "fastqc is not installed, and is needed for performing sample quality control"

def test_fastq_screen_installed():
    assert application_installed('fastq_screen'), "fastq_screen is not installed, and is needed for performing sample quality control"

