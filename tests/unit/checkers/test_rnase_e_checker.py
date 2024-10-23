import pytest
from genedesign.checkers.rnase_e_checker import RNaseEChecker

@pytest.fixture
def rnase_e_checker():
    checker = RNaseEChecker()
    checker.initiate()
    return checker

def test_no_rnase_e_sites(rnase_e_checker):
    # TODO: IMPLEMENT
    pass
