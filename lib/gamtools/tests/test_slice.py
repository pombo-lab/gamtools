import gamtools.slice_wrapper

def test_slice_pkg_structure():
    assert hasattr(gamtools.slice_wrapper, "slice")

def test_slice_compilation():
    assert gamtools.slice_wrapper.test_compilation() == 4

