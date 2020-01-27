#!/bin/env python3
import pytest
import sys
from os import getcwd, path, pardir
sys.path.append(path.abspath(path.join(getcwd())))
import smored

def test_bad_fasta_path():
	__reads__ = False
	__unclassified__ = False
	smored.read_processor("tests", 35, "no_such_sample", None, None)

def test_bad_fasta_content():
	__reads__ = False
	__unclassified__ = False
	smored.read_processor("tests", 35, "test_sample", None, None)
