#!/bin/env python3
import pytest
import sys
from os import getcwd, path, pardir
sys.path.append(path.abspath(path.join(getcwd(), pardir)))
from abil_urdocaller import reverse_complement
def test_reverse_complement_01():
    assert reverse_complement('ATCG') == 'CGAT'
def test_reverse_complement_02():
    assert reverse_complement('ATCG') != 'TAGC'
def test_reverse_complement_03():
    assert reverse_complement('TTAA') != 'AATT'

