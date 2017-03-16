import io
import pytest
from numpy.testing import assert_array_equal
import numpy as np
from gamtools import call_windows, segregation

fixture_two_samples = io.StringIO(
u"""chrom   start   stop    Sample_A    Sample_B
chr1    0   50000   4   0
chr1    50000   100000  0   5
chr1    100000  150000  3   0
chr1    150000  200000  0   0
chr1    200000  250000  0   6
""")

data_two_samples = segregation.open_segregation(fixture_two_samples)

def test_fixed_threshold_4():

    threshold_function = call_windows.fixed_threshold_fitting_func(4)

    fitting_result = threshold_function(data_two_samples.Sample_A)

    assert fitting_result['read_threshold'] == 4
    assert 'counts' in fitting_result
    assert 'breaks' in fitting_result
    assert fitting_result['params'] is None


def test_fixed_coverage_thresholding():

    threshold_function = call_windows.fixed_threshold_fitting_func(4)

    segregation_table, fitting_result = call_windows.do_coverage_thresholding(
        data_two_samples, None, threshold_function)

    assert_array_equal(segregation_table.Sample_A,
                       np.array([0,0,0,0,0]))

    assert_array_equal(segregation_table.Sample_B,
                       np.array([0,1,0,0,1]))

