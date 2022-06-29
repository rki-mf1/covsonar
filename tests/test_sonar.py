import pytest

from covsonar import sonar


def test_setup_db(tmp_path):
    parsed_args = sonar.parse_args(['setup', '--db', str(tmp_path / 'test.db')])
    retval = sonar.main(parsed_args)
    assert retval == 0
