import re
from genedesign.models.transcript import Transcript
from genedesign.models.rbs_option import RBSOption
from genedesign.models.operon import Operon

class RNaseEChecker:
    """
    Description: 
    This class checks sequences for RNase E cleavage sites, which generally come in the form
    (A/G)AUU(A/U)
    
    Input (run method):
        seq (str): An operon DNA sequence

    Output:
        Tuple[bool, int, List[int]]: A tuple containing:
            - rnase_e_site_present (bool): False if an RNase E cleavage site is present
            - num_sites (int): Number of cleavage sites present in the sequence
            - positions (List[int]): List of indices where the cleavage sites are present
    """

    def __init__(self) -> None:
        self.pattern = None

    def initiate(self) -> None:
        """
        Defines the regex pattern used to check the input sequence.
        """

        self.pattern = re.compile(r"[AU]AUU[AU]")

    def run(self, seq):
        """
        Checks the input sequence for cleavage sites based on the RNase E pattern.

        Parameters:
            seq (str): An operon DNA sequence.
        
        Returns: 
            Tuple containing:
                - rnase_e_site_present (bool)
                - num_sites(int)
                - positions (List[int])
                - matches (List[str])
            As defined in the class definition. 
        """

        matches = list(self.pattern.finditer(seq))

        positions = [match.start() for match in matches]

        num_sites = len(positions)

        if num_sites > 0:
            rnase_e_site_present = False
        else:
            rnase_e_site_present = True

        return rnase_e_site_present, num_sites, positions, matches
    
if __name__ == "__main__":
    """
    Main method for running the RNaseEChecker on a hardcoded Operon.
    This method initializes RNaseEChecker, runs it on a hardcoded Operon, 
    and prints the results including a boolean indicating if an
    RNase E site is present, the number of sites if they exist, and the position(s) 
    of those sites. 
    """

    rnase_e_checker = RNaseEChecker()
    rnase_e_checker.initiate()

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
    operon = Operon(
        transcripts=[transcript1, transcript2],
        promoter="TTGACA",
        terminator="TGA"
    )

    operon_seq = operon.promoter + "".join([transcript.rbs.utr + transcript.rbs.cds for transcript in operon.transcripts]) + operon.terminator

    site_in_operon, num_operon_sites, operon_site_pos, operon_matches = rnase_e_checker.run(operon_seq)

    print(f"No RNase E cleavage sites are present in the operon: {site_in_operon}")
    print(f"Number of sites in operon: {num_operon_sites}")
    print(f"The operon sites are located at: {operon_site_pos}")
    print(f"Matching sequences within the operon: {operon_matches}")
