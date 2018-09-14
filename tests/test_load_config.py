#!/bin/env python3
import pytest
import sys
from os import getcwd, path, pardir
sys.path.append(path.abspath(path.join(getcwd())))
import abil_urdocaller
def test_good_config():
	abil_urdocaller.load_config('tests/config_good.txt')
	assert True != False
	assert abil_urdocaller.__config_dict__ == {'loci': {'alleles': 'tests/repDB-rename.fasta'}, 'profile': {'profile': 'tests/profile.txt'}}
