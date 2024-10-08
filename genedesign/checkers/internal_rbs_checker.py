import re 
from genedesign.models.transcript import Transcript
from genedesign.models.rbs_option import RBSOption
from genedesign.models.operon import Operon

class InternalRBSChecker:
    """
    Description: This class checks sequences for internal RBS binding sites
    containing a Shine-Dalgarno sequence and a start codon 7-9 bases away that
    could lead to unintended translation initiation. 

    Inputs (run method):
        seq (str): An operon DNA sequence

    Output: 
        Tuple[bool, int, List[int]]: a tuple containing:
            - internal_rbs_present (bool): False if an internal RBS is present.
            - num_sites (int): Number of RBSes present in the sequence.
            - positions (List[int]): List of indices where the RBSes are present.
    """

    def __init__(self) -> None:
        """
        Defines the regex pattern used to check the input for a Shine-Dalgarno 
        sequence and a start codon.
        """

        self.rbs_pattern = None

    def initiate(self) -> None:
        """
        Defines the regex pattern used to check the input sequence.
        """
        
        self.rbs_pattern = re.compile(r"AUUAUU[AUCG]{7,9}AUG")

    def run(self, seq):
        """
        Checks the input sequence for an RBS based on the pattern.

        Parameters:
            seq (str): An operon DNA sequence

        Returns: 
            Tuple containing:
                - internal_rbs_present (bool)
                - num_sites (int)
                - positions (List[int])
            As defined in the class definition. 
        """

        seq.replace('T', 'U')

        matches = list(self.rbs_pattern.finditer(seq))

        positions = [match.start() for match in matches]

        num_sites = len(positions)

        if num_sites > 0:
            internal_rbs_present = False
        else:
            internal_rbs_present = True

        return internal_rbs_present, num_sites, positions
    
if __name__ == "__main__":
    """
    Main method for running the InternalRBSChecker on a hardcoded Operon. 
    This method initializes InternalRBSChecker, runs it on a hardcoded Operon,
    annd prints the results including a boolean indicating if an internal RBS is
    present, the number of sites if they exist, and the position(s) of those sites.
    """

    internal_rbs_checker = InternalRBSChecker()
    internal_rbs_checker.initiate()

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
    operon = Operon(
        transcripts=[transcript1, transcript2],
        promoter="TTGACA",
        terminator="TGA"
    )

    operon_seq = operon.promoter + "".join([transcript.rbs.utr + transcript.rbs.cds for transcript in operon.transcripts]) + operon.terminator

    internal_rbs, num_sites, positions = internal_rbs_checker.run(operon_seq)

    print(f"No internal RBSes are present in the operon: {internal_rbs}")
    print(f"Number of sites in operon: {num_sites}")
    print(f"The operon sites are located at: {positions}")