def exclude_ignored(self, ignores):
        """
        Excludes RBSOptions that are in the ignores set.

        Parameters: 
            ignores (Set[[RBSOptions]): Set of RBSOptions to be ignored.

        Returns:
            List[RBSOption]: List of RBSOptions that are not in the ignores set.
        """     

        return [option for option in self.rbs_options if option not in ignores]