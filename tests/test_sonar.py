import pytest

from covsonar import sonar


def test_help():
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parsed_args = sonar.parse_args(['--help'])
        sonar.main(parsed_args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_setup_db(tmp_path):
    parsed_args = sonar.parse_args(['setup', '--db', str(tmp_path / 'test.db')])
    retval = sonar.main(parsed_args)
    assert retval == 0
