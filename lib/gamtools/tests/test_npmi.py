from gamtools import segregation, cosegregation
import io
from numpy.testing import assert_array_equal, assert_array_almost_equal
import pytest
import numpy as np



fixture_invalid_data = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K
chr3    0       50000   0 0 0 0 1 2 1 0 0 0 0
chr3    50000   100000  0 0 0 0 0 1 1 0 1 0 0
""")

data_invalid_data = segregation.open_segregation(fixture_invalid_data)

fixture_region_c = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K
chr3    0       50000   0 0 0 0 1 1 1 0 0 0 0
chr3    50000   100000  0 0 0 0 0 1 1 0 1 0 0
""")

data_region_c = segregation.open_segregation(fixture_region_c)

fixture_region_d = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0
""")

data_region_d = segregation.open_segregation(fixture_region_d)

fixture_region_e = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   0 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0
""")

data_region_e = segregation.open_segregation(fixture_region_e)

fixture_region_f = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0
""")

data_region_f = segregation.open_segregation(fixture_region_f)

fixture_region_g = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1 1 1 1
""")

data_region_g = segregation.open_segregation(fixture_region_g)

fixture_region_h = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
""")

data_region_h = segregation.open_segregation(fixture_region_h)

fixture_region_i = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1
""")

data_region_i = segregation.open_segregation(fixture_region_i)

fixture_region_k = io.StringIO(
u"""chrom start  stop    A B C D E F G H I J K L M N O P Q R S T
chr3    0       50000   0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
""")

data_region_k = segregation.open_segregation(fixture_region_k)

def test_npmi_positive():

    npmi_res = cosegregation.get_npmi_from_regions(data_region_d, data_region_e)

    assert_array_almost_equal(npmi_res, np.array([[ 0.5129416 ]]))

def test_npmi_zero():

    npmi_res = cosegregation.get_npmi_from_regions(data_region_d, data_region_f)

    assert_array_almost_equal(npmi_res, np.array([[ 0.0 ]]))

def test_npmi_negative():

    npmi_res = cosegregation.get_npmi_from_regions(data_region_d, data_region_g)

    assert_array_almost_equal(npmi_res, np.array([[ -0.39794 ]]))

def test_npmi_max():

    npmi_res = cosegregation.get_npmi_from_regions(data_region_d, data_region_d)

    assert_array_almost_equal(npmi_res, np.array([[ 1.0 ]]))

#def test_npmi_min():
#
    #npmi_res = cosegregation.get_npmi_from_regions(data_region_i, data_region_k)

    #assert_array_almost_equal(npmi_res, np.array([[ -1.0 ]]))

def test_npmi_not_detected():

    npmi_res = cosegregation.get_npmi_from_regions(data_region_d, data_region_h)

    assert_array_almost_equal(npmi_res, np.array([[ np.nan ]]))

def test_npmi_multiple_windows():

    npmi_res = cosegregation.get_npmi_from_regions(data_region_c)

    assert_array_almost_equal(npmi_res, np.array([[ 1.0, 0.5243108 ],
                                                     [ 0.5243108, 1.0 ]]))

def test_npmi_invalid_data():

   with pytest.raises(cosegregation.InvalidDataError):
       cosegregation.get_npmi_from_regions(data_invalid_data)
