try:
    from StringIO import StringIO
except:
    from io import StringIO
import unittest

import pytest

from gamtools.qc.merge import merge_stats


@pytest.fixture
def stats_file_1():
    return StringIO(
u"""Sample	col1	col2
001	1	11
002	2	12
""")

@pytest.fixture
def stats_file_2():
    return StringIO(
u"""Sample	col3	col4
001	21	31
002	22	32
""")

@pytest.fixture
def stats_file_3():
    return StringIO(
u"""Sample	col3	col4
a01	21	31
b02	22	32
""")

@pytest.fixture
def stats_file_4():
    return StringIO(
u"""Sample	col3	col4
a01	21	31
b02	22	32
""")

def test_string_sample_names(stats_file_3, stats_file_4):

    output_buffer = StringIO()
    merge_stats([stats_file_3, stats_file_4], output_buffer)

    output_buffer.seek(0)

    first_line = [row.split('\t')[0] for row in output_buffer.readlines()]

    assert first_line[1] == 'a01'
    assert first_line[2] == 'b02'

def test_numeric_sample_names(stats_file_1, stats_file_2):

    output_buffer = StringIO()
    merge_stats([stats_file_1, stats_file_2], output_buffer)

    output_buffer.seek(0)

    first_line = [row.split('\t')[0] for row in output_buffer.readlines()]

    assert first_line[1] == '001'
    assert first_line[2] == '002'
