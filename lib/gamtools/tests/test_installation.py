import os

import pytest
from gamtools import main

def test_tests_installed():
    assert os.path.exists(main.test_directory)

def test_examples_installed():
    assert os.path.exists(main.get_example('qc_parameters.example.cfg'))

def test_scripts_installed():
    assert os.path.exists(main.get_script('mapping_stats.sh'))
