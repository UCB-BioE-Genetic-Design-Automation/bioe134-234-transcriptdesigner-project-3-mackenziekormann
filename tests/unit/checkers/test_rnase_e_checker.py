import pytest
from genedesign.models.transcript import Transcript
from genedesign.models.rbs_option import RBSOption
from genedesign.models.operon import Operon
from genedesign.checkers.rnase_e_checker import RNaseEChecker

@pytest.fixture
def checker():
    checker = RNaseEChecker()
    checker.initiate()
    return checker

@pytest.fixture
def sample_operon():
    rbs1 = RBSOption(
        utr="AGGGAUUAAGGAUAUUAU", 
        cds="AUGGCUAAGGAGGAUGCUAA", 
        gene_name="gene1", 
        first_six_aas="MAKGMC"
    )
    rbs2 = RBSOption(
        utr="AUGGAUUAGGAUUUAA", 
        cds="AUGGUAUUUAAGGCUAGG", 
        gene_name="gene2", 
        first_six_aas="MGYAGG"
    )
    transcript1 = Transcript(
        rbs=rbs1,
        peptide="MAGMC",
        codons=["AUG", "GCU", "AAG", "GAG", "GCU"]
    )
    transcript2 = Transcript(
        rbs=rbs2,
        peptide="MGYAG",
        codons=["AUG", "GUA", "UUU", "AAG", "GCU"]
    )
    return Operon(
        transcripts=[transcript1, transcript2],
        promoter="TTGACA",
        terminator="TGA"
    )

def test_no_rnase_e_site(checker):
    # A sequence with no RNase E site
    seq_no_rnase_e = "TTGACAGGCCGGAAGGCUUCGAAUGCGAUGCGAGCGGAACGCGCGA"
    rnase_e_site_present, num_sites, positions, matches = checker.run(seq_no_rnase_e)

    # Expecting no RNase E sites to be found
    assert rnase_e_site_present == True  # True means no RNase E site found
    assert num_sites == 0
    assert positions == []
    assert matches == []

def test_rnase_e_site_present(checker):
    # Sequence with a single RNase E site
    seq_with_rnase_e = "AUAUAUUAUUGCGCUAUUGCAUG"
    rnase_e_site_present, num_sites, positions, matches = checker.run(seq_with_rnase_e)

    # Expecting one RNase E site to be found
    assert rnase_e_site_present == False
    assert num_sites == 1
    assert positions == [3]  # Expected start position of RNase E site
    assert matches[0].group() == "AUUAU"

def test_multiple_rnase_e_sites(checker):
    # Sequence with multiple RNase E sites
    seq_with_multiple_rnase_e = "GGAUUAUGAUUAUUUAUAUGCAUAUAUUAUGCGCUU"
    rnase_e_site_present, num_sites, positions, matches = checker.run(seq_with_multiple_rnase_e)

    # Expecting multiple RNase E sites to be found
    assert rnase_e_site_present == False
    assert num_sites == 3
    assert positions == [3, 10, 24]  # Expected positions of RNase E sites
    assert matches[0].group() == "AUUAU"
    assert matches[1].group() == "AUUAU"
    assert matches[2].group() == "AUUAU"

def test_empty_sequence(checker):
    # Empty sequence test
    empty_seq = ""
    rnase_e_site_present, num_sites, positions, matches = checker.run(empty_seq)

    # Expecting no RNase E sites to be found
    assert rnase_e_site_present == True
    assert num_sites == 0
    assert positions == []
    assert matches == []

def test_boundary_case_no_match(checker):
    # Sequence close to matching but without the exact RNase E pattern
    boundary_seq = "AAUUAU"
    rnase_e_site_present, num_sites, positions, matches = checker.run(boundary_seq)

    # Expecting no RNase E sites
    assert rnase_e_site_present == True
    assert num_sites == 0
    assert positions == []
    assert matches == []

def test_sample_operon_sequence(checker, sample_operon):
    # Run RNaseEChecker on a constructed operon
    operon_seq = sample_operon.promoter + "".join(
        [transcript.rbs.utr + transcript.rbs.cds for transcript in sample_operon.transcripts]
    ) + sample_operon.terminator

    rnase_e_site_present, num_sites, positions, matches = checker.run(operon_seq)

    # Check if RNase E sites are correctly identified
    assert rnase_e_site_present == False  # There should be RNase E sites
    assert num_sites > 0  # Ensure at least one RNase E site exists
    assert len(positions) == num_sites
    assert all(isinstance(pos, int) for pos in positions)

