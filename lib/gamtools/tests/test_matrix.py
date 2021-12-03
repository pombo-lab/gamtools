import io

import pytest
import numpy as np
from gamtools import matrix


@pytest.fixture
def small_matrix():
    return io.StringIO(
u"""	chr1:0-1000	chr1:1000-2000	chr1:2000-3000	chr1:3000-4000
chr1:0-1000	1	2	3	4
chr1:1000-2000	5	6	7	8
chr1:2000-3000	9	10	11	12
chr1:3000-4000	13	14	15	16
""")


@pytest.fixture
def trans_matrix():
    return io.StringIO(
u"""	chr2:0-1000	chr2:1000-2000	chr2:2000-3000	chr2:3000-4000
chr1:0-1000	1	2	3	4
chr1:1000-2000	5	6	7	8
chr1:2000-3000	9	10	11	12
chr1:3000-4000	13	14	15	16
""")


def test_matrix_object_subregion(small_matrix):

    matrix_obj = matrix.read_txt(small_matrix)
    loc_string = 'chr1:2000-4000'
    (w1, w2), subregion = matrix.region_from_locations(matrix_obj, loc_string)

    np.testing.assert_array_equal(subregion,
                       np.array([[11,12],[15,16]]))


def test_matrix_object_offcentre_region(small_matrix):

    matrix_obj = matrix.read_txt(small_matrix)
    loc_string1 = 'chr1:2000-4000'
    loc_string2 = 'chr1:0-2000'

    (w1, w2), subregion = matrix.region_from_locations(matrix_obj, loc_string1, loc_string2)

    np.testing.assert_array_equal(subregion,
                       np.array([[9,10],[13,14]]))


def test_trans_matrix_object_offcentre_region(trans_matrix):

    matrix_obj = matrix.read_txt(trans_matrix)
    loc_string1 = 'chr1:2000-4000'
    loc_string2 = 'chr2:0-2000'

    (w1, w2), subregion = matrix.region_from_locations(matrix_obj, loc_string1, loc_string2)

    np.testing.assert_array_equal(subregion,
                       np.array([[9,10],[13,14]]))


def test_matrix_file_subregion(small_matrix):

    loc_string = 'chr1:2000-4000'
    (w1, w2), subregion = matrix.open_region_from_locations(
        small_matrix, loc_string, file_type='txt')

    np.testing.assert_array_equal(subregion,
                       np.array([[11,12],[15,16]]))


def test_matrix_file_offcentre_region(small_matrix):

    loc_string1 = 'chr1:2000-4000'
    loc_string2 = 'chr1:0-2000'

    (w1, w2), subregion = matrix.open_region_from_locations(
        small_matrix, loc_string1, loc_string2, file_type='txt')

    np.testing.assert_array_equal(subregion,
                       np.array([[9,10],[13,14]]))


def test_trans_matrix_file_offcentre_region(trans_matrix):

    loc_string1 = 'chr1:2000-4000'
    loc_string2 = 'chr2:0-2000'

    (w1, w2), subregion = matrix.open_region_from_locations(
        trans_matrix, loc_string1, loc_string2, file_type='txt')

    np.testing.assert_array_equal(subregion,
                       np.array([[9,10],[13,14]]))
