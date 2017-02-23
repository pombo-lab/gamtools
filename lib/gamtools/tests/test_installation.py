import os

import pytest
from gamtools import main

def test_tests_installed():
    assert os.path.exists(main.test_directory)
