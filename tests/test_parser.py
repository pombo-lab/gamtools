import pytest
from gamtools import main


def test_no_args():
    with pytest.raises(SystemExit) as error_msg:
        print(main.parser.parse_args([]))

