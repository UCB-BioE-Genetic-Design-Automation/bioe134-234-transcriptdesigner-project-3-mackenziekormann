from seq_utils.calc_edit_distance import calculate_edit_distance

def compare_peptides(self, rbs_option, cds):
        """
        Compares the peptide sequence of the input CDS and the RBSOption source gene.

        Parameters: 
            rbs_option (RBSOption): The RBSOption to compare to the given CDS.
            cds (str): The CDS to compare to the given RBSOption. 

        Returns:
            int: The edit distance between the peptide sequences of the RBS and the CDS.
        """

        cds_peptide = self.translator.run(cds[:18])
        return calculate_edit_distance(cds_peptide, rbs_option.first_six_aas)