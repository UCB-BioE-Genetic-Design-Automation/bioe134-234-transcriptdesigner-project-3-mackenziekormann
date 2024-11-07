import pytest
from genedesign.models.transcript import Transcript
from genedesign.models.rbs_option import RBSOption
from genedesign.models.operon import Operon
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

@pytest.fixture
def checker():
    checker = InternalRBSChecker()
    checker.initiate()
    return checker

@pytest.fixture
def sample_operon():
    # Sample operon sequences as described in the initial code
    rbs1 = RBSOption(
        utr="AUUAUUACGGAUG",  
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

def test_no_internal_rbs(checker, sample_operon):
    operon_seq = sample_operon.promoter + "".join(
        [transcript.rbs.utr + transcript.rbs.cds for transcript in sample_operon.transcripts]
    ) + sample_operon.terminator
    internal_rbs, num_sites, positions, matches = checker.run(operon_seq)

    # Assertions to verify there is no internal RBS
    assert internal_rbs == True  # True means no internal RBS found
    assert num_sites == 0
    assert positions == []
    assert matches == []

def test_internal_rbs_present(checker):
    # Sequence with an internal RBS to trigger detection
    seq_with_internal_rbs = "TTGACAAUUAUUAAAGGUAUGGGCTGA"
    internal_rbs, num_sites, positions, matches = checker.run(seq_with_internal_rbs)

    # Assertions to verify the presence of internal RBS
    assert internal_rbs == False  # False indicates an internal RBS was found
    assert num_sites == 1
    assert positions == [6]  # Position where the RBS starts
    assert matches[0].group() == "AUUAUUAAAGGUAUG"

def test_multiple_internal_rbs_sites(checker):
    # Sequence with multiple internal RBS sites
    seq_with_multiple_rbs = "TTGACAAUUAUUAAAGGUAUGGGCAUUAUUCCGGAUGTGACGA"
    internal_rbs, num_sites, positions, matches = checker.run(seq_with_multiple_rbs)

    # Assertions for multiple RBS sites
    assert internal_rbs == False
    assert num_sites == 2
    assert positions == [6, 20]  # Expected starting positions of RBS patterns
    assert matches[0].group() == "AUUAUUAAAGGUAUG"
    assert matches[1].group() == "AUUAUUCCGGAUG"

def test_empty_sequence(checker):
    # Test with an empty sequence
    empty_seq = ""
    internal_rbs, num_sites, positions, matches = checker.run(empty_seq)

    # Assertions for empty input
    assert internal_rbs == True
    assert num_sites == 0
    assert positions == []
    assert matches == []

def test_short_sequence_no_match(checker):
    # Test with a sequence that is too short to contain an RBS
    short_seq = "AUUA"
    internal_rbs, num_sites, positions, matches = checker.run(short_seq)

    # Assertions for short sequence
    assert internal_rbs == True
    assert num_sites == 0
    assert positions == []
    assert matches == []

def test_boundary_conditions(checker):
    # Test with a sequence just at the boundary of matching criteria
    boundary_seq = "AUUAUUAAAGG"  # Pattern just below expected distance for RBS
    internal_rbs, num_sites, positions, matches = checker.run(boundary_seq)

    # Assertions for boundary condition (should not match as it lacks the full RBS)
    assert internal_rbs == True
    assert num_sites == 0
    assert positions == []
    assert matches == []

