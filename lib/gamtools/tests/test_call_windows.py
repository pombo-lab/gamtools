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

fixture_cov_two_samples = io.StringIO(
u"""chrom   start   stop    Sample_A    Sample_B
chr1    0   50000   0.002   0
chr1    50000   100000  0   0.004
chr1    100000  150000  0.002   0
chr1    150000  200000  0   0
chr1    200000  250000  0   0.002
""")

fixture_cov_two_samples_two_reads = io.StringIO(
u"""chrom   start   stop    Sample_A    Sample_B
chr1    0   50000   0.004   0
chr1    50000   100000  0   0.002
chr1    100000  150000  0.004   0
chr1    150000  200000  0   0
chr1    200000  250000  0   0.004
""")

fixture_no_positive_windows = np.array(
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

fixture_all_positive_windows = np.array(
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

fixture_first_window = np.array(
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

fixture_last_window = np.array(
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])

fixture_no_orphans = np.array(
    [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0])

fixture_1_orphan = np.array(
    [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0])

fixture_all_orphans = np.array(
    [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0])

fixture_coverage_data = np.array(
    [0, 0, 0, 0, 60, 60, 60, 50, 60, 60, 60, 60, 50, 60, 50, 0, 0, 0, 0, 0,
     0, 45, 0, 40, 0, 35, 0, 30, 0, 25, 0, 20, 0, 15, 0, 10, 0, 5, 0, 0])

data_two_samples = segregation.open_segregation(fixture_two_samples)

data_cov_two_samples = segregation.open_segregation(fixture_cov_two_samples)

data_cov_two_samples_two_reads = segregation.open_segregation(fixture_cov_two_samples_two_reads)

def test_fixed_threshold_4():

    threshold_function = call_windows.fixed_threshold_fitting_func(4)

    fitting_result = threshold_function(data_two_samples.Sample_A)

    assert fitting_result['coverage_threshold'] == 4.0
    assert fitting_result['params'] is None


def test_fixed_coverage_thresholding():

    threshold_function = call_windows.fixed_threshold_fitting_func(4)

    segregation_table, fitting_result = call_windows.do_coverage_thresholding(
        data_two_samples, None, threshold_function)

    assert_array_equal(segregation_table.Sample_A,
                       np.array([0,0,0,0,0]))

    assert_array_equal(segregation_table.Sample_B,
                       np.array([0,1,0,0,1]))


def test_orphan_windows():

    assert call_windows.pct_orphan_windows(fixture_no_positive_windows) == np.Inf

    assert call_windows.pct_orphan_windows(fixture_all_positive_windows) == 0

    assert call_windows.pct_orphan_windows(fixture_first_window) == 1.

    assert call_windows.pct_orphan_windows(fixture_last_window) == 1.

    assert call_windows.pct_orphan_windows(fixture_no_orphans) == 0.

    assert call_windows.pct_orphan_windows(fixture_1_orphan) == (1./ 10.)

    assert call_windows.pct_orphan_windows(fixture_all_orphans) == 1.


def test_orphan_threshold_no_clip():

    window_func = call_windows.orphan_windows_fitting_func('1,99', clip=False)

    results = window_func(fixture_coverage_data)

    assert results['coverage_threshold'] > 45
    assert results['coverage_threshold'] < 50


def test_orphan_threshold_clipping():

    window_func = call_windows.orphan_windows_fitting_func('1,99', clip=True)

    results = window_func(fixture_coverage_data)

    assert results['coverage_threshold'] == 60


def test_infer_read_length():

    inferred_length = call_windows.infer_single_read_coverage(data_cov_two_samples)

    assert inferred_length == 0.002

    inferred_length = call_windows.infer_single_read_coverage(data_cov_two_samples_two_reads)

    assert inferred_length == 0.002
