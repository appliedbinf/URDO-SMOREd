#!/bin/env python3
import pytest
import sys
from os import getcwd, path, pardir
#sys.path.append(path.abspath("~/repo"))
sys.path.append(path.abspath(path.join(getcwd())))
print(sys.path)
from abil_urdocaller import get_fasta_dict
def test_get_fasta_dict_singleline():
	assert get_fasta_dict("tests//test_get_fasta_dict.fasta") == {'seq1': {'sequence': 'ATGC'}}

def test_get_fasta_dict_multiline():
	assert get_fasta_dict("tests/test_get_fasta_dict_multiline.fasta") ==	{'multilineseq': {'sequence': 'ATCGTGCA'}}

def test_get_fasta_dict_bad_fasta(): 
	assert get_fasta_dict("tests/test_get_fasta_dict_bad.fasta") != {'seq1': {'sequence': 'ATGC'}, 'multilineseq': {'sequence': 'ATCGTGCA'}}
	assert get_fasta_dict("tests/test_get_fasta_dict_bad.fasta") == {'seq1': {'sequence': 'ATGCmultilineseqATCGTGCA'}}

