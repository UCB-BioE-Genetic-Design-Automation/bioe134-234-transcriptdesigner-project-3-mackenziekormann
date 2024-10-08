from models.rbs_option import RBSOption
from Bio import SeqIO
from collections import defaultdict
from seq_utils.Translate import Translate
from seq_utils.reverse_complement import reverse_complement
from seq_utils.exclude_ignored import exclude_ignored
from seq_utils.check_secondary_structure import check_secondary_structure
from seq_utils.compare_peptides import compare_peptides

class RBSChooser:
    """
    A simple RBS selection algorithm that chooses an RBS from a list of options, excluding any RBS in the ignore set.
    """

    def __init__(self):
        self.rbs_options = []
        self.translator = Translate()
        self.translator.initiate()
        self.gene_data = None
        self.pruned_data = None
        self.merged_data = None

    @staticmethod
    def extract_genes_info(genbank_file):
        """
        Extracts gene information from a given genbank file. 

        Parameters:
            genbank_file: Path to the file containing sequence information.

        Returns:
            gene_dict (dict): Dictionary mapping locus ID to gene name, UTR, CDS< and sequence.
        """

        gene_dict = defaultdict(dict)
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                if feature.type == "gene":
                    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                    gene_name = feature.qualifiers.get("gene", [None])[0]

                    cds_feature = None
                    for cds in record.features:
                        if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                            cds_feature = cds
                            break

                    if cds_feature:
                        start, end = cds_feature.location.start, cds_feature.location.end
                        strand = cds_feature.location.strand
                        if strand == 1: 
                            utr_start = max(0, start - 50)
                            utr_seq = record.seq[utr_start:start]
                        else:
                            utr_start = end
                            utr_seq = reverse_complement(record.seq[utr_start:utr_start + 50])

                        cds_seq = cds_feature.extract(record.seq)
                        gene_dict[locus_tag] = {
                            "gene": gene_name,
                            "UTR": utr_seq,
                            "CDS": cds_seq
                        }
        return gene_dict
    
    @staticmethod
    def prune_gene_info(proteomics_file):
        """
        Prunes proteomics data to only maintain top 5% of entries by abundance. 

        Parameters: 
            proteomics_file: Path to the file containing proteomics data.

        Returns: 
            pruned_data (dict): A dictionary containing the pruned proteomics data. 
        """

        proteomics_data = {}

        with open(proteomics_file, 'r') as file:
            for line in file:
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.strip().split('\t')

                if len(parts) != 2:
                    continue

                protein_id = parts[0].split('.')[1]
                abundance = float(parts[1])

                proteomics_data[protein_id] = abundance

        sorted_proteomics_data = dict(sorted(proteomics_data.items(), key=lambda item: item[1], reverse=True))
        top_5_percent = max(1, int(0.05 * len(sorted_proteomics_data)))
        pruned_data = dict(list(sorted_proteomics_data.items())[:top_5_percent])

        return pruned_data
    
    @staticmethod
    def merge_data(pruned_data, gene_info):
        """
        Creates a dictionary where the key is the locus taga nd the value is the RBS/CDS data based on proteomics data. 

        Parameters: 
            pruned_data (dict): The pruned proteomics data containing IDs and abundances.
            gene_info (dict): The gene information dictionary containing UTR, RBS, CDS, and gene names. 

        Returns:
            merged_dict (dict): A dictionary where the key is the locus tag and the value is a dictionary with RBS and CDS.
        """  

        merged_data = {}

        for id in pruned_data.keys():
            if id in gene_info.keys():
                rbs = gene_info[id].get("UTR")
                cds = gene_info[id].get("CDS")

                merged_data[id] = {
                    "RBS": rbs,
                    "CDS": cds
                }    
        
        return merged_data

    def initiate(self) -> None:
        """
        Initializes RBSChooser with RBS options. 
        """

        self.gene_data = self.extract_genes_info('genedesign/data/sequence.gb')
        self.pruned_data = self.prune_gene_info('genedesign/data/511145-WHOLE_ORGANISM-integrated.txt')
        self.merged_data = self.merge_data(self.pruned_data, self.gene_data)

        for locus_tag, gene_info in self.gene_data.items():
            if locus_tag in self.pruned_data:
                utr = gene_info.get('UTR')
                cds = gene_info.get('CDS')
                gene_name = gene_info.get('gene')
                first_six_aas = self.translator.run(str(cds)[:18])

                rbs_option = RBSOption(utr, cds, gene_name, first_six_aas)
                self.rbs_options.append(rbs_option)
        
    def run(self, cds: str, ignores: set) -> RBSOption:
        """
        Selects an RBS that is not in the ignore set, accounting for secondary structures and sequence similarity.
        
        Parameters:
            cds (str): The coding sequence.
            ignores (set): A set of RBS options to ignore.
        
        Returns:
            RBSOption: The selected RBS option.
        """

        filtered_rbs_options = exclude_ignored(self, ignores)

        best_option = None
        best_score = float('inf')

        for option in filtered_rbs_options:
            hairpins = check_secondary_structure(self, option, cds)[0]
            distance = compare_peptides(self, option, cds)
            score = hairpins + distance
            
            if score < best_score:
                best_option = option
                best_score = score

        return best_option

if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT"

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)
