import pytest
import numpy as np
from gamtools import main


def test_no_args():
    with pytest.raises(SystemExit) as error_msg:
        print(main.parser.parse_args([]))


def test_fixed_threshold():

    args = main.parser.parse_args(['call_windows', '-x', '4', '/dev/null'])

    fitting_result = args.fitting_function(np.arange(100))

    assert fitting_result['coverage_threshold'] == 4
    assert fitting_result['params'] is None


def test_fixed_threshold_5():

    args = main.parser.parse_args(['call_windows', '-x', '5', '/dev/null'])

    fitting_result = args.fitting_function(np.arange(100))

    assert fitting_result['coverage_threshold'] == 5.0
    assert fitting_result['params'] is None


