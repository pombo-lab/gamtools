from GamTools import segmentation
import StringIO
from nose.tools import assert_raises
from numpy.testing import assert_array_equal
import numpy as np

fixture_two_samples = StringIO.StringIO("""chrom   start   stop    Sample_A    Sample_B
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

data_two_samples = segmentation.open_segmentation(fixture_two_samples)

fixture_two_windows = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J
chr1    0       50000   0 0 0 1 1 1 1 1 1 1
chr1    50000   100000  0 1 1 0 0 0 1 1 1 1
""")

data_two_windows = segmentation.open_segmentation(fixture_two_windows)

fixture_three_windows = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e f g h i j
chr1    0       50000   0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
chr1    50000   100000  0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
chr1    100000  150000  0 1 1 0 0 0 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1
""")

data_three_windows = segmentation.open_segmentation(fixture_three_windows)

fixture_invalid_data = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J
chr1    0       50000   0 0 0 2 1 1 1 1 1 1
chr1    50000   100000  0 1 1 0 0 0 1 1 1 1
""")

data_invalid_data = segmentation.open_segmentation(fixture_invalid_data)

#########################################
#
# segmentation.index_from_interval tests
#
#########################################

def test_interval_within_bin():
    interval = 'chr1', 50100, 50300
    start_index, stop_index = segmentation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    assert len(found_windows) == 1
    print found_windows
    found_chrom, found_start, found_stop = found_windows[0]
    assert found_chrom == 'chr1'
    assert found_start == 50000
    assert found_stop == 100000

def test_interval_is_bin():
    interval = 'chr1', 50000, 100000
    start_index, stop_index = segmentation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    assert len(found_windows) == 1
    print found_windows
    found_chrom, found_start, found_stop = found_windows[0]
    assert found_chrom == 'chr1'
    assert found_start == 50000
    assert found_stop == 100000

def test_interval_overlaps_two_bins():
    interval = 'chr1', 50500, 100500
    start_index, stop_index = segmentation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    print found_windows
    assert len(found_windows) == 2
    assert found_windows[0] == ('chr1', 50000, 100000)
    assert found_windows[-1] == ('chr1', 100000, 150000)

def test_interval_overlaps_many_bins():
    interval = 'chr1', 50500, 300500
    start_index, stop_index = segmentation.index_from_interval(data_two_samples,
                                                               interval)
    found_windows = data_two_samples.index[start_index:stop_index]
    print found_windows
    assert len(found_windows) == 6
    assert found_windows[0] == ('chr1', 50000, 100000)
    assert found_windows[-1] == ('chr1', 300000, 350000)

def test_interval_end_before_start():
    interval = 'chr1', 300500, 50500 
    assert_raises(ValueError,
                  segmentation.index_from_interval,
                  data_two_samples,
                  interval)

def test_invalid_chromosome():
    interval = 'chr3', 50000, 100000
    assert_raises(segmentation.InvalidChromError,
                  segmentation.index_from_interval,
                  data_two_samples,
                  interval)

#########################################
#
# segmentation.index_from_interval tests
#
#########################################

def test_cosegregation_two_loci():
    
    segregation_freqs = segmentation.cosegregation_frequency(np.array(data_two_windows))

    assert_array_equal(segregation_freqs, np.array([[ 1., 2.],
                                                    [ 3., 4.]]))

def test_cosegregation_three_loci():
    
    segregation_freqs = segmentation.cosegregation_frequency(np.array(data_three_windows))

    assert_array_equal(segregation_freqs, np.array([[[ 1., 2.],
                                                     [ 3., 4.]],
                                                    [[ 5., 6.],
                                                     [ 7., 8.]]]))

def test_cosegregation_invalid_data():

    assert_raises(IndexError,
                  segmentation.cosegregation_frequency,
                  np.array(data_invalid_data))

def test_index_combinations_one_region():

    regions = (['a'] * 3,)

    combinations = segmentation.get_index_combinations(regions)

    print list(combinations)

    assert len(list(combinations)) == 9

def test_index_combinations_three_regions():

    regions = (['a'], ['b'] * 2, ['c'] * 3)

    combinations = segmentation.get_index_combinations(regions)

    assert len(list(combinations)) == 6
