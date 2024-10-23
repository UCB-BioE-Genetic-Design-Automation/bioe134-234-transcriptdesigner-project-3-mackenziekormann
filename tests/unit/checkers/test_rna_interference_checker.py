import pytest
from genedesign.checkers.rna_interference_checker import RNAInterferenceChecker

@pytest.fixture
def rna_interference_checker():
    checker = RNAInterferenceChecker()
    checker.initiate()
    return checker

def test_no_interference(rna_interference_checker):
    # TODO: IMPLEMENT
    pass

