from collections import namedtuple
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

import gamtools.slice_wrapper
import gamtools.gam_slice
import gamtools.main

Namespace = namedtuple('Namespace', (
                            'segregation_file_path',
                            'output_dir',
                            'skip_chroms',
                            'slice_thickness',
                            'nuclear_radius',
                            'genome_size'))

def test_slice_pkg_structure():
    assert hasattr(gamtools.slice_wrapper, "slice")

def test_slice_compilation():
    assert gamtools.slice_wrapper.test_compilation() == 4

@patch('gamtools.gam_slice.get_slice_output_dirs', return_value=['a', 'b', 'c', 'd', 'e'])
@patch('gamtools.gam_slice.slice_wrapper.slice')
def test_slice_no_G(slice_mock, output_dirs_mock):

    segregation_path = gamtools.main.get_example('example_segregation.table')
    args = Namespace(segregation_path, '/path/to/matrix',
                     ['chrY'], 0.22, 100, None)

    gamtools.gam_slice.run_slice_from_args(args)

    print(slice_mock)

    slice_mock.assert_called_with('e', 'a', 'b', 'c', 'd', 5, 24000, 1000, 0.22, 100)


@patch('gamtools.gam_slice.get_slice_output_dirs', return_value=['a', 'b', 'c', 'd', 'e'])
@patch('gamtools.gam_slice.slice_wrapper.slice')
def test_slice_with_G(slice_mock, output_dirs_mock):

    segregation_path = gamtools.main.get_example('example_segregation.table')
    args = Namespace(segregation_path, '/path/to/matrix',
                     ['chrY'], 0.22, 100, 1e9)

    gamtools.gam_slice.run_slice_from_args(args)

    slice_mock.assert_called_with('e', 'a', 'b', 'c', 'd', 5, 1e9, 1000, 0.22, 100)
