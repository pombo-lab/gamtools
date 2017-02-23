from gamtools import segregation
import io
import pytest

fixture_two_samples = io.StringIO(
u"""chrom   start   stop    Sample_A    Sample_B
chr1    0   50000   0   0
chr1    50000   100000  0   0
chr1    100000  150000  0   0
chr1    150000  200000  0   0
chr1    200000  250000  0   0
chr1    250000  300000  0   0
chr1    300000  350000  0   0
chr1    350000  400000  0   0
chr1    400000  450000  0   0
""")

data_two_samples = segregation.open_segregation(fixture_two_samples)

#########################################
#
# segregation.index_from_interval tests
#
#########################################

def test_interval_within_bin():
    interval = 'chr1', 50100, 50300
    start_index, stop_index = segregation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    assert len(found_windows) == 1
    print(found_windows)
    found_chrom, found_start, found_stop = found_windows[0]
    assert found_chrom == 'chr1'
    assert found_start == 50000
    assert found_stop == 100000

def test_interval_is_bin():
    interval = 'chr1', 50000, 100000
    start_index, stop_index = segregation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    assert len(found_windows) == 1
    print(found_windows)
    found_chrom, found_start, found_stop = found_windows[0]
    assert found_chrom == 'chr1'
    assert found_start == 50000
    assert found_stop == 100000

def test_interval_overlaps_two_bins():
    interval = 'chr1', 50500, 100500
    start_index, stop_index = segregation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    print(found_windows)
    assert len(found_windows) == 2
    assert found_windows[0] == ('chr1', 50000, 100000)
    assert found_windows[-1] == ('chr1', 100000, 150000)

def test_interval_overlaps_many_bins():
    interval = 'chr1', 50500, 300500
    start_index, stop_index = segregation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    print(found_windows)
    assert len(found_windows) == 6
    assert found_windows[0] == ('chr1', 50000, 100000)
    assert found_windows[-1] == ('chr1', 300000, 350000)

def test_interval_end_before_start():
    interval = 'chr1', 300500, 50500 
    with pytest.raises(ValueError):
        segregation.index_from_interval(data_two_samples, interval)

def test_invalid_chromosome():
    interval = 'chr3', 50000, 100000
    with pytest.raises(segregation.InvalidChromError):
        segregation.index_from_interval(data_two_samples, interval)

