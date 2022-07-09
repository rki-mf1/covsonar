import os
from sonar import sonarBasics
import pytest
import re

def test_setup_but_file_exists(get_tmpfile_name):
	fname = get_tmpfile_name
	open(fname, "w").close()
	with pytest.raises(SystemExit) as pytest_wrapped_e:
		sonarBasics.setup_db(fname)
	assert pytest_wrapped_e.type == SystemExit
	assert re.match(r'^setup error: .* does already exist.$', pytest_wrapped_e.value.code) != None
