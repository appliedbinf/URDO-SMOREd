#!/bin/env python3
import pytest
import sys
from os import getcwd, path, pardir
sys.path.append(path.abspath(path.join(getcwd())))
import smored
def test_good_config():
	smored.load_config('tests/config_good.txt')
	assert True != False
	assert smored.__config_dict__ == {'loci': {'alleles': 'tests/repDB-rename.fasta'}, 'profile': {'profile': 'tests/profile.txt'}}

# def test_bad_config():
# 	with pytest.raises(SystemExit) as pytest_wrapped_e:
# 		abil_urdocaller.load_config('tests/config_bad.txt')
# 	print(pytest_wrapped_e)
# 	assert pytest_wrapped_e.type == SystemExit
# 	assert pytest_wrapped_e.value.code == 1

