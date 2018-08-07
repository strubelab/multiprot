"""
Script that contains the higher level implementation of the ranch wrapper,
to handle single and multiple chain scenarios

"""


##### FOR THE EARLY IMPLEMENTATION OF MULTIPROT, WE WILL ASSUME ONLY ONE CHAIN
##### IN THE INPUT


class Builder:
    """
    Class that will contain all the necessary methods to call the Ranch and
    pulchra wrappers multiple times and build single and multiple chain 
    structures. This class will be used by the multiprot script.
    
    """

    def __init__(self, chains):
        """
        :param args: Object that contains the arguments parsed from the command line
        :type args: argparse.Namespace object created by calling parser.parse_args()
        """
        self.chains = chains

    def build(self):
        """
        Method that will build the models and return a return a list with 10
        PDBModels
        """

        for i in range(len(chains)):
            if chains[i][4] == False:   # If it is still not modeled
                call = r.Ranch(*chains[i][0], chains=chains[i][1], symmetry=args.symmetry, 
                    fixed = chains[i][2], symtemplate = chains[i][3], 
                    pool_sym=args.pool_sym)
                models = call.run()
                chains[i][4] = True
                
                ## QUESTION: TAKE ONE OR TEN MODELS FOR NEXT STEP?

            # Search every ('file.pdb', 'chainID') pair in other chains of chain[x][5] to 
            # see if the same multichain pdb in chain i is used in another chain j
            for pdb, chainID in chains[i][5]:
                # For each pair of ('file.pdb', 'chainID') in present chain
                
                for j in range(i+1, len(chains)):
                    # For each of the remaining chains in 'chains'
                    
                    for _pdb, _chainID in chains[j][5]:
                        # For each ('file.pdb', 'chainID') pair in chain j[5]
                        
                        if pdb == _pdb:
                            # The same pdb is used, probably with a different chain
                            # Now what
                            # I'll tell  you what
                            # Extract the chain with new coordinates, embed chains[i],
                            # EXCEPT parts that are in chains[j] bound to chains[i],
                            # those will be taken from chains[i] and passed to ranch to
                            # be modeled as individual fixed domains in chains[j]






