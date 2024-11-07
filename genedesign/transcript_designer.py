import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.transcript_to_seq import transcript_to_seq
from genedesign.checkers.codon_checker import CodonChecker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence, chooses an RBS, 
    and optimizes codon selection using sliding window and guided random approaches.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.codonChecker = CodonChecker()
        self.forbiddenChecker = ForbiddenSequenceChecker()

    def initiate(self) -> None:
        """
        Initializes the codon table, RBS chooser, codon checker, and forbidden sequence checker.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        self.codonChecker.initiate()
        self.forbiddenChecker.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': ["GCG"], 'C': ["TGC"], 'D': ["GAT"], 'E': ["GAA"], 'F': ["TTC"],
            'G': ["GGT"], 'H': ["CAC"], 'I': ["ATC"], 'K': ["AAA"], 'L': ["CTG"],
            'M': ["ATG"], 'N': ["AAC"], 'P': ["CCG"], 'Q': ["CAG"], 'R': ["CGT"],
            'S': ["TCT"], 'T': ["ACC"], 'V': ["GTT"], 'W': ["TGG"], 'Y': ["TAC"]
        }

    def run(self, peptide: str, ignores: set) -> Transcript:
        codons = []
        window_size = 3  # Sliding window size for amino acids

        for i in range(0, len(peptide), window_size):
            # Extract a segment (window) of the peptide sequence
            window = peptide[i:i + window_size]
            
            # Generate candidate solutions for the current window
            candidate_solutions = [
                [random.choice(self.aminoAcidToCodon[aa]) for aa in window] for _ in range(10)
            ]

            # Filter candidate solutions by forbidden sequences
            valid_solutions = []
            for solution in candidate_solutions:
                codons_in_solution = ''.join(solution)
                
                # Check with forbidden sequence checker
                if self.forbiddenChecker.run(codons_in_solution)[0]:  # True means no forbidden sequences found
                    valid_solutions.append(solution)

            # Select the first valid solution or fallback if none
            print(valid_solutions[0])
            if valid_solutions:
                best_solution = valid_solutions[0]
            else:
                # Fallback: Use the first generated solution if no valid ones
                best_solution = candidate_solutions[0]
            
            codons.extend(best_solution)

        # Add a stop codon to finalize the CDS
        codons.append("TAA")
        cds = ''.join(codons)

        # Choose RBS using ignores set
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the finalized transcript
        transcript = Transcript(selectedRBS, peptide, codons)
        return transcript


if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)

