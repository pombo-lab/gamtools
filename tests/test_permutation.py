try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from gamtools import segregation, permutation
import numpy as np
import io

fixture_two_samples = io.StringIO(
u"""chrom   start   stop    Sample_A    Sample_B
chr1    0       50000   1   0
chr1    50000   100000  1   0
chr1    100000  150000  0   1
chr1    150000  200000  0   0
chr2    200000  250000  1   1
chr2    250000  300000  0   1
chr2    300000  350000  1   0
chr2    350000  400000  0   0
chr2    400000  450000  0   1
""")

data_two_samples = segregation.open_segregation(fixture_two_samples)


def test_permute_by_1():

    sample_a_perm_1 = permutation.permute_by_offset(data_two_samples.Sample_A, 1)
    expected_a = np.array([0, 1, 1, 0, 0, 1, 0, 1, 0])

    print('obs:', sample_a_perm_1)
    print('exp:', expected_a)

    assert np.array_equal(sample_a_perm_1,
                          expected_a)


def test_permute_by_2():

    sample_a_perm_2 = permutation.permute_by_offset(data_two_samples.Sample_A, 2)
    expected_a = np.array([0, 0, 1, 1, 0, 0, 1, 0, 1])

    print('obs:', sample_a_perm_2)
    print('exp:', expected_a)

    assert np.array_equal(sample_a_perm_2,
                          expected_a)


def test_permute_by_more_than_length():

    sample_a_perm_2 = permutation.permute_by_offset(data_two_samples.Sample_A, 11)
    expected_a = np.array([0, 0, 1, 1, 0, 0, 1, 0, 1])

    print('obs:', sample_a_perm_2)
    print('exp:', expected_a)

    assert np.array_equal(sample_a_perm_2,
                          expected_a)


def test_permute_by_multiple_of_length():

    sample_a_perm_2 = permutation.permute_by_offset(data_two_samples.Sample_A, 29)
    expected_a = np.array([0, 0, 1, 1, 0, 0, 1, 0, 1])

    print('obs:', sample_a_perm_2)
    print('exp:', expected_a)

    assert np.array_equal(sample_a_perm_2,
                          expected_a)


@patch('gamtools.permutation.np.random.randint')
def test_only_permute_mappable(mock_randint):

    mock_randint.return_value = 1

    segregation_perm_1 = permutation.permute_segregation(data_two_samples)
    expected_a = np.array([0, 1, 1, 0, 0, 1, 0, 0, 1])

    print(segregation_perm_1.Sample_A.values)
    print(expected_a)

    assert np.array_equal(segregation_perm_1.Sample_A.values,
                          expected_a)


@patch('gamtools.permutation.np.random.randint')
def test_each_np_permuted_separately(mock_randint):

    mock_randint.return_value = 1

    segregation_perm_1 = permutation.permute_segregation(data_two_samples)
    expected_a = np.array([0, 1, 1, 0, 0, 1, 0, 0, 1])
    expected_b = np.array([1, 0, 0, 0, 1, 1, 1, 0, 0])

    print('obs A:', segregation_perm_1.Sample_A.values)
    print('exp A:', expected_a)
    print('obs B:', segregation_perm_1.Sample_B.values)
    print('exp B:', expected_b)

    assert np.array_equal(segregation_perm_1.Sample_A.values,
                          expected_a)

    assert np.array_equal(segregation_perm_1.Sample_B.values,
                          expected_b)


