#!/bin/env python3
import pytest
import sys
from os import getcwd, path, pardir
sys.path.append(path.abspath("~/repo"))
from abil_urdocaller import reverse_complement
def test_reverse_complement():
    assert reverse_complement('ATCG') == 'CGAT'
    assert reverse_complement('aNt') == 'TNA'
    assert reverse_complement('ATCG') != 'TAGC'
    assert reverse_complement('TTAA') != 'AATT'
