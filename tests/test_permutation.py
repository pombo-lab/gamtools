from gamtools import segregation, permutation
import io
import pytest

fixture_two_samples = io.StringIO(
u"""chrom   start   stop    Sample_A    Sample_B
chr1    0   50000   1   1
chr1    50000   100000  1   1
chr1    100000  150000  0   0
chr1    150000  200000  0   0
chr2    200000  250000  1   1
chr2    250000  300000  0   0
chr2    300000  350000  1   1
chr2    350000  400000  0   0
chr2    400000  450000  0   0
""")

data_two_samples = segregation.open_segregation(fixture_two_samples)


