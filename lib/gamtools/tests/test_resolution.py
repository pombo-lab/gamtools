from gamtools import segregation, resolution, utils
import io
import pytest

fixture_table = io.StringIO(
u"""chrom   start   stop    Sample_1    Sample_2    Sample_3    Sample_4    Sample_5    Sample_6    Sample_7    Sample_8    Sample_9    Sample_10    Sample_11    Sample_12    Sample_13    Sample_14    Sample_15    Sample_16    Sample_17    Sample_18    Sample_19    Sample_20    Sample_21    Sample_22
chr1    0      8000         0           1           0           0           0           0           0           0           0           0            0            0            0            0            0            0            0            0            0            0            0            0
chr1    8000   16000        0           0           0           0           0           0           0           0           0           0            0            0            0            0            0            0            0            0            0            0            0            0
""")

data_table = segregation.open_segregation(fixture_table)

bigger_table = segregation.open_segregation(
    utils.get_example('example_segregation.table'))

def test_detection_efficiency_visible():
    eff = resolution.get_detection_efficiency(data_table, slice_thickness = 200, nuclear_radius = 10000, genome_size = 1e9)
    assert 0.78937754 < eff < 0.78937755

def test_detection_efficiency_all():
    eff = resolution.get_detection_efficiency(data_table, slice_thickness = 200, nuclear_radius = 10000, genome_size = 1e9, only_visible=False)
    assert 0.39239383 < eff < 0.39239384

def test_window_coverage():
    cov = resolution.get_segregation_coverage_qc(bigger_table)
    assert cov == 2

def test_window_coverage_quantile():
    cov = resolution.get_segregation_coverage_qc(bigger_table, quantile=0.6)
    assert cov == 3

def test_window_coverage_invisible():
    cov = resolution.get_segregation_coverage_qc(bigger_table, only_visible=False)
    assert cov == 0
