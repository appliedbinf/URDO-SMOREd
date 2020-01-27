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

def test_bad_config_path():
	with pytest.raises(OSError) as pytest_wrapped_e:
		smored.load_config('tests/config_bad.txt')
	print(pytest_wrapped_e)
	assert pytest_wrapped_e.type == OSError

def se():
	raise SystemExit(1)

# @pytest.mark.skip
def test_bad_config_noprofile():
	smored.load_config('tests/config_noprofile.txt')

