from seq_utils.hairpin_counter import hairpin_counter

def check_secondary_structure(self, rbs_option, cds):
        """
        Checks for the potential secondary structure formation in the combined RBS-CDS sequence.

        Parameters: 
            rbs_option (RBSOption): The RBSOption to check with the given CDS.
            cds (str): The CDS sequence to check with the given RBSOption.

        Returns:
            int: The number of secondary structures (hairpins) detected in the RBS + CDS.
        """

        combined_seq = rbs_option.utr + cds
        return hairpin_counter(combined_seq)