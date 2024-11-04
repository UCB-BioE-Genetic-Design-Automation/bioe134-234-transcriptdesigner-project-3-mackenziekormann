import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.hairpin_checker import hairpin_checker
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

    def initiate(self) -> None:
        """
        Initializes the codon table, RBS chooser, and codon checker.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        self.codonChecker.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': ["GCG"], 'C': ["TGC"], 'D': ["GAT"], 'E': ["GAA"], 'F': ["TTC"],
            'G': ["GGT"], 'H': ["CAC"], 'I': ["ATC"], 'K': ["AAA"], 'L': ["CTG"],
            'M': ["ATG"], 'N': ["AAC"], 'P': ["CCG"], 'Q': ["CAG"], 'R': ["CGT"],
            'S': ["TCT"], 'T': ["ACC"], 'V': ["GTT"], 'W': ["TGG"], 'Y': ["TAC"]
        }

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA, selects an RBS, and optimizes codon selection.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with selected RBS and optimized codons.
        """
        codons = []
        window_size = 3  # 3 amino acids per window

        for i in range(0, len(peptide), window_size):
            window = peptide[i:i + window_size]

            # Generate candidate solutions for each window using guided random selection
            candidate_solutions = [
                [random.choice(self.aminoAcidToCodon[aa]) for aa in window] for _ in range(10)
            ]

            # Filter and rank candidate solutions based on CodonChecker results
            ranked_solutions = []
            for solution in candidate_solutions:
                codons_in_solution = ''.join(solution)
                check_result = self.codonChecker.run(solution)
                
                if check_result[0]:  # Passes the thresholds
                    ranked_solutions.append((check_result, solution))

            # Choose the best solution if any valid candidates exist; else use fallback
            if ranked_solutions:
                best_solution = max(ranked_solutions, key=lambda x: (x[0][1], -x[0][2], x[0][3]))  # High diversity, low rare codons, high CAI
                codons.extend(best_solution[1])
            else:
                # Fallback to highest CAI codons for each amino acid in the window
                fallback_solution = [self.aminoAcidToCodon[aa][0] for aa in window]
                codons.extend(fallback_solution)

        # Append the stop codon (TAA in this case)
        codons.append("TAA")

        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Check the transcript for secondary structures and forbidden sequences
        transcript = Transcript(selectedRBS, peptide, codons)
        if not hairpin_checker(transcript_to_seq(transcript).upper())[0]:
            raise ValueError("Hairpin(s) found in the sequence.")

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

