"""
Script that contains the higher level implementation of the ranch wrapper

"""

#RG: I think this should also be a class with methods, defaults and all
#RG: the actual script should then be very small and only deal with command line parsing and argument cleanup
#JG:  okay so I will divide this in two then? One only for argument parsing and
#     other for the methods

import argparse
import biskit as b
import random
import ranch as r

# If type=divide does not work, try action='append_const'
#RG: no big deal but instead of `string`, which looks like a type definition, convention is `s`
def divide(s):  
    return tuple(s.split(':'))
    # if the string is not properly formatted
    # raise argparse.ArgumentTypeError(msg)

# Supported symmetries for --symmetry argument
supp_sym = ['p'+str(i) for i in range(1,20)]
supp_sym += ['p'+str(i) for i in range(22,132,10)] + ['p222']

# Add positional and optional arguments
# NOTES: try formatter_class=argparse.MetavarTypeHelpFormatter
#        if no argument is given, show help (and possibly raise error)
#        look at customizing file parsing
#        look at exiting methods

parsero = argparse.ArgumentParser(usage='', 
    description='''Build multiple chain-proteins. The arguments can also
    be read from a file, in which case the file name must have the @ prefix''')

parsero.add_argument('--chain', '-c', action='append', nargs='+', type=divide, 
    help='Add a new chain to the model', required=True)

parsero.add_argument('--split', '-spl', default=None, 
    help='Split one or more chains from a given PDB')

parsero.add_argument('--symmetry', '-sym', default='p1',
    help='What kind of symmetry do you wish to have in your molecule. Supported\n\
    symmetries are: p1, p2, …, p19 (nineteen-fold), p22, p32, p42, p52, p62,\n\
    …, p122, p222.', choices=supp_sym)

parsero.add_argument('--symtemplate', '-t', default=[], action='append', 
    help='Which domain will be the symmetry core, in case of symmetry other than\
     p1 specified')

## Check the options, because ranch wrapper only supports full words right now
parsero.add_argument('--pool_sym', '-o', default='mix', 
    choices=['mixed', 'm', 'symmetry', 's', 'asymmetry', 'a'], 
    help='Specify the overall symmetry of the molecules to be produced, i.e. \
    all symmetric [s], all asymmetric [a] or mixed. [m]')

parsero.add_argument('--fixed', '-f', default=[], nargs='*',
    help='Specify one or more domains to be fixed in their original coordinates.')

parsero.add_argument('args', nargs=argparse.REMAINDER)

#argument_default = argparse.SUPPRESS

args = parsero.parse_args()
print(args)

# parsero.Namespace() ?

# vars(args)  returns dictionary with attributes

## Create PDBModels

chains = []    # Arguments for ranch for each chain
multidom = []  # List of tuples [('1234.pdb', 'A'), ... ] for each chain
                    # for all the multichain domains with chains specified

class Multiprot:
    """
    Class that will contain all the necessary methods to call the Ranch wrapper
    multiple times and build multiple chain structures. This class will only be
    used by the XXXX script.
    
    """

    chains = []    # Arguments for ranch for each chain

    def __init__(self, args):
        """
        :param args: Object that contains the arguments parsed from the command line
        :type args: argparse.Namespace object created by calling parser.parse_args()
        """
        self.args = args

    def create_args(self):
        """
        Method that takes the arguments parsed from the command line and 
        """
        for i in range(len(self.args.chain)):    # For each chain
            rdomains = []
            rchains = {}
            rfixed = []
            rsymtemp = None
            rmulti = []

            rdomains = [b.PDBModel(dom) for dom in self.args.chain[i][j][0] \
                    for j in range(len(self.args.chain[i]))]
            rchains = {b.PDBModel(dom):chain for dom, chain in 
                    self.args.chain[i][j][0], self.args.chain[i][j][1] \
                    if len(self.args.chain[i][j])==2 \
                    for j in range(len(self.args.chain[i]))}
            rfixed = [b.PDBModel(dom) for dom in self.args.chain[i][j][0] \
                    if self.args.chain[i][j][0] in self.args.fixed \
                    for j in range(len(self.args.chain[i]))]
            rsymtemp = [b.PDBModel(dom) for dom in self.args.chain[i][j][0] \
                    ]

            chains.append([[b.PDBmodel(dom) for dom in self.args.chain[i][j][0]], 
                {b.PDBModel(dom):chain for dom, chain in self.args.chain[i][j][0],
                 self.args.chain[i][j][1]}])

            for j in range(len(self.args.chain[i])):    # For each component of the chain

                if self.args.chain[i][j][0][-4:]=='.pdb':
                    pdb = b.PDBModel(self.args.chain[i][j][0])
                        if len(self.args.chain[i][j])==2:  # If chain to be taken is specified
                            rchains[pdb] = self.args.chain[i][j][1]
                            rmulti.append((self.args.chain[i][j][0], self.args.chain[i][j][1]))
                        if self.args.chain[i][j][0] in self.args.fixed:  # If domain will be fixed
                            rfixed.append(pdb)
                        if self.args.chain[i][j][0] in self.args.symtemplate:  # If is symtemplate
                            rsymtemp = pdb    # THERE CAN ONLY BE ONE

                    rdomains.append(pdb)
                    continue
                rdomains.append(self.args.chain[i][j][0])

            # chains[i] = (domains/linkers, chains dict, symtemplate, already modeled)
            # rmulti is a list of tuples [('1234.pdb', 'A'), ... ] for each domain with
            # a chain specified
            chains.append((rdomains, rchains, rfixed, rsymtemp, False, rmulti))

        random.shuffle(chains)

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






