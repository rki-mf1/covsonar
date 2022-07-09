import re

import pytest

from covsonar.basics import sonarBasics


def test_setup_but_file_exists(tmpfile_name):
    fname = tmpfile_name
    open(fname, "w").close()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sonarBasics.setup_db(fname)
    assert pytest_wrapped_e.type == SystemExit
    assert (
        re.match(r"^setup error: .* does already exist.$", pytest_wrapped_e.value.code)
        is not None
    )
