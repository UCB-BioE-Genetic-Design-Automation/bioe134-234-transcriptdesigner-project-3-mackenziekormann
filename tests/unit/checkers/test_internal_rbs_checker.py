import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

@pytest.fixture
def internal_rbs_checker():
    checker = InternalRBSChecker()
    checker.initiate()
    return checker

def test_no_internal_rbs(internal_rbs_checker):
    # TODO: IMPLEMENT
    pass

def test_one_internal_rbs(internal_rbs_checker):
    # TODO: IMPLEMENT
    pass

def test_multiple_internal_rbs(internal_rbs_checker):
    # TODO: IMPLEMENT
    pass

def test_edge_case(internal_rbs_checker):
    # TODO: IMPLEMENT
    pass