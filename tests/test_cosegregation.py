from GamTools import segmentation, cosegregation
import StringIO
from numpy.testing import assert_array_equal, assert_raises, assert_array_almost_equal
import numpy as np
from mock import patch

fixture_window1_only = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J
chr1    0       50000   0 0 0 1 1 1 1 1 1 1
chr1    50000   100000  0 0 0 0 0 0 0 0 0 0
""")

data_window1_only = segmentation.open_segmentation(fixture_window1_only)

fixture_window2_only = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J
chr1    0       50000   0 0 0 0 0 0 0 0 0 0
chr1    50000   100000  0 1 1 0 0 0 1 1 1 1
""")

data_window2_only = segmentation.open_segmentation(fixture_window2_only)

fixture_region_a = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K
chr1    0       50000   1 1 1 1 1 1 0 0 0 0 0
chr1    50000   100000  0 0 0 1 1 1 1 1 0 0 0
chr1    100000  150000  0 0 0 0 0 1 1 0 1 1 0
""")

data_region_a = segmentation.open_segmentation(fixture_region_a)

fixture_region_b = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K
chr2    0       50000   0 0 0 0 1 1 1 0 0 0 0
chr2    50000   100000  0 0 0 0 0 1 1 0 1 0 0
chr2    100000  150000  0 0 0 0 0 0 0 1 1 0 1
""")

data_region_b = segmentation.open_segmentation(fixture_region_b)

fixture_region_c = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K
chr3    0       50000   0 0 0 0 1 1 1 0 0 0 0
chr3    50000   100000  0 0 0 0 0 1 1 0 1 0 0
""")

data_region_c = segmentation.open_segmentation(fixture_region_c)

fixture_invalid_data = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K
chr3    0       50000   0 0 0 0 1 2 1 0 0 0 0
chr3    50000   100000  0 0 0 0 0 1 1 0 1 0 0
""")

data_invalid_data = segmentation.open_segmentation(fixture_invalid_data)

#########################################
#
# segmentation.cosegregation_frequency tests
#
#########################################

def test_cosegregation_one_region():

    segregation_freqs = cosegregation.get_cosegregation_from_regions(data_region_a)

    assert_array_equal(segregation_freqs, np.array([[ 6., 3., 1.],
                                                    [ 3., 5., 2.],
                                                    [ 1., 2., 4.]]))

def test_cosegregation_two_regions():

    segregation_freqs = cosegregation.get_cosegregation_from_regions(data_region_a,
                                                                     data_region_b)

    assert_array_equal(segregation_freqs, np.array([[ 2., 1., 0.],
                                                    [ 3., 2., 1.],
                                                    [ 2., 3., 1.]]))

def test_cosegregation_three_regions():

    segregation_freqs = cosegregation.get_cosegregation_from_regions(data_region_c,
                                                                     data_region_c,
                                                                     data_region_c)

    assert_array_equal(segregation_freqs, np.array([[[ 3., 2. ],
                                                     [ 2., 2. ]],
                                                    [[ 2., 2. ],
                                                     [ 2., 3. ]]]))

def test_cosegregation_missing_windows():

    segregation_freqs = cosegregation.get_cosegregation_from_regions(data_window1_only,
                                                                     data_window2_only)

    assert_array_equal(segregation_freqs, np.array([[ 0., 4.],
                                                    [ 0., 0.]]))

def test_cosegregation_invalid_data():

    assert_raises(cosegregation.InvalidDataError,
                  cosegregation.get_cosegregation_from_regions,
                  data_invalid_data)

def test_cosegregation_min():

    cosegregation_res = cosegregation.get_cosegregation_from_regions(data_region_i, data_region_k)

    assert_array_equal(cosegregation_res, np.array([[ 0.0 ]]))


#########################################
#
# segmentation.linkage tests
#
#########################################

fixture_region_d = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0
""")

data_region_d = segmentation.open_segmentation(fixture_region_d)

fixture_region_e = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   0 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0
""")

data_region_e = segmentation.open_segmentation(fixture_region_e)

fixture_region_f = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0
""")

data_region_f = segmentation.open_segmentation(fixture_region_f)

fixture_region_g = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1 1 1 1
""")

data_region_g = segmentation.open_segmentation(fixture_region_g)

fixture_region_h = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
""")

data_region_h = segmentation.open_segmentation(fixture_region_h)

fixture_region_i = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1
""")

data_region_i = segmentation.open_segmentation(fixture_region_i)

fixture_region_j = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0
""")

data_region_j = segmentation.open_segmentation(fixture_region_j)

fixture_region_k = StringIO.StringIO(
"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
""")

data_region_k = segmentation.open_segmentation(fixture_region_k)

def test_linkage_positive():

    linkage_res = cosegregation.get_linkage_from_regions(data_region_d, data_region_e)

    assert_array_almost_equal(linkage_res, np.array([[ 0.15 ]]))

def test_linkage_zero():

    linkage_res = cosegregation.get_linkage_from_regions(data_region_d, data_region_f)

    assert_array_almost_equal(linkage_res, np.array([[ 0.0 ]]))

def test_linkage_negative():

    linkage_res = cosegregation.get_linkage_from_regions(data_region_d, data_region_g)

    assert_array_almost_equal(linkage_res, np.array([[ -0.15 ]]))

def test_linkage_min():

    linkage_res = cosegregation.get_linkage_from_regions(data_region_i, data_region_k)

    assert_array_almost_equal(linkage_res, np.array([[ -0.25 ]]))

# Make sure we don't fail the test because of the warning message
@patch('GamTools.cosegregation.warnings.warn')
def test_3d_linkage_positive(mock_warnings):

    linkage_res = cosegregation.get_linkage_from_regions(data_region_d,
                                                         data_region_f,
                                                         data_region_i)

    assert_array_almost_equal(linkage_res, np.array([[[ 0.125 ]]]))

# Make sure we don't fail the test because of the warning message
@patch('GamTools.cosegregation.warnings.warn')
def test_3d_linkage_zero(mock_warnings):

    linkage_res = cosegregation.get_linkage_from_regions(data_region_d,
                                                         data_region_e,
                                                         data_region_f)

    assert_array_almost_equal(linkage_res, np.array([[[ 0.0 ]]]))

# Make sure we don't fail the test because of the warning message
@patch('GamTools.cosegregation.warnings.warn')
def test_3d_linkage_negative(mock_warnings):

    linkage_res = cosegregation.get_linkage_from_regions(data_region_f,
                                                         data_region_j,
                                                         data_region_k)

    assert_array_almost_equal(linkage_res, np.array([[[ -0.1 ]]]))

@patch('GamTools.cosegregation.warnings.warn')
def test_3d_linkage_warning(mock_warnings):

    linkage_res = cosegregation.get_linkage_from_regions(data_region_f,
                                                         data_region_j,
                                                         data_region_k)

    assert mock_warnings.called

def test_linkage_not_detected():

    linkage_res = cosegregation.get_linkage_from_regions(data_region_d, data_region_h)

    assert_array_almost_equal(linkage_res, np.array([[ np.nan ]]))

def test_linkage_multiple_windows():

    linkage_res = cosegregation.get_linkage_from_regions(data_region_c)

    assert_array_almost_equal(linkage_res, np.array([[ 0.198347, 0.107438 ],
                                                     [ 0.107438, 0.198347 ]]))

def test_linkage_invalid_data():

    assert_raises(cosegregation.InvalidDataError,
                  cosegregation.get_linkage_from_regions,
                  data_invalid_data)



#########################################
#
# segmentation.linkage tests
#
#########################################

def test_dprime_positive():

    dprime_res = cosegregation.get_dprime_from_regions(data_region_d, data_region_e)

    assert_array_almost_equal(dprime_res, np.array([[ 0.6 ]]))

def test_dprime_zero():

    dprime_res = cosegregation.get_dprime_from_regions(data_region_d, data_region_f)

    assert_array_almost_equal(dprime_res, np.array([[ 0.0 ]]))

def test_dprime_negative():

    dprime_res = cosegregation.get_dprime_from_regions(data_region_d, data_region_g)

    assert_array_almost_equal(dprime_res, np.array([[ -0.6 ]]))

def test_dprime_max():

    dprime_res = cosegregation.get_dprime_from_regions(data_region_d, data_region_d)

    assert_array_almost_equal(dprime_res, np.array([[ 1.0 ]]))

def test_dprime_min():

    dprime_res = cosegregation.get_dprime_from_regions(data_region_i, data_region_k)

    assert_array_almost_equal(dprime_res, np.array([[ -1.0 ]]))

def test_dprime_not_detected():

    dprime_res = cosegregation.get_dprime_from_regions(data_region_d, data_region_h)

    assert_array_almost_equal(dprime_res, np.array([[ np.nan ]]))

def test_dprime_multiple_windows():

    dprime_res = cosegregation.get_dprime_from_regions(data_region_c)

    assert_array_almost_equal(dprime_res, np.array([[ 1.0, 0.541667 ],
                                                     [ 0.541667, 1.0 ]]))

def test_dprime_invalid_data():

    assert_raises(cosegregation.InvalidDataError,
                  cosegregation.get_dprime_from_regions,
                  data_invalid_data)
