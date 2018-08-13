"""
Script that handles command line argument parsing and calls Builder module

"""

## Question: prevent the user to 

#RG: I think this should also be a class with methods, defaults and all
#RG: the actual script should then be very small and only deal with command line parsing and argument cleanup
#JG:  okay so I will divide this in two then? One only for argument parsing and
#     other for the methods

import argparse
import biskit as b
import random
import builder
import os

# If type=divide does not work, try action='append_const'
#RG: no big deal but instead of `string`, which looks like a type definition, convention is `s`
def divide(s):  
    return tuple(s.split(':'))
    # if the string is not properly formatted
    # raise argparse.ArgumentTypeError(msg)

def path_exists(s):
    if not os.path.exists(s):
        raise argparse.ArgumentTypeError('Invalid destination for output files.')
    return s

# Supported symmetries for --symmetry argument
supp_sym = ['p'+str(i) for i in range(1,20)]
supp_sym += ['p'+str(i) for i in range(22,132,10)] + ['p222']

# Add positional and optional arguments
# NOTES: try formatter_class=argparse.MetavarTypeHelpFormatter
#        if no argument is given, show help (and possibly raise error)
#        look at customizing file parsing
#        look at exiting methods

parser = argparse.ArgumentParser(usage='', 
    description='''Build multiple chain-proteins. The arguments can also
    be read from a file, in which case the file name must have the @ prefix''')

parser.add_argument('--chain', '-c', action='append', nargs='+', type=divide, 
    help='Add a new chain to the model', required=True)

parser.add_argument('--split', '-spl', default=None, 
    help='Split one or more chains from a given PDB')

parser.add_argument('--symmetry', '-sym', default='p1',
    help='What kind of symmetry do you wish to have in your molecule. Supported\n\
    symmetries are: p1, p2, …, p19 (nineteen-fold), p22, p32, p42, p52, p62,\n\
    …, p122, p222.', choices=supp_sym)

parser.add_argument('--symtemplate', '-t', default=[], action='append', 
    help='Which domain will be the symmetry core, in case of symmetry other than\
     p1 specified')

## Check the options, because ranch wrapper only supports full words right now
parser.add_argument('--poolsym', '-o', default='mix', 
    choices=['mixed', 'm', 'symmetry', 's', 'asymmetry', 'a'], 
    help='Specify the overall symmetry of the molecules to be produced, i.e. \
    all symmetric [s], all asymmetric [a] or mixed. [m]')

parser.add_argument('--fixed', '-f', default=[], nargs='*',
    help='Specify one or more domains to be fixed in their original coordinates.')

parser.add_argument('--destination', '-d', default=os.getcwd(), type=path_exists, 
    help='Specify the directory where the output models will be saved (default cwd)')

parser.add_argument('--debug', action='store_true')

# parser.add_argument('args', nargs=argparse.REMAINDER, help="Additional key=value\
#      parameters are passed on to 'Executor.__init__'. For example:\n\
#         debug   - 0|1, keep all temporary files (default: 0)")

#argument_default = argparse.SUPPRESS

args = parser.parse_args()

# vars(args)  returns dictionary with attributes


def create_args(args):
    """
    Function that takes the arguments parsed from the command line and
    returns a list with

    :param args:    Namespace object created by calling .parse_args() method on 
                    the command line arguments  
    """
    
    chains = []

    for i in range(len(args.chain)):    # For each chain
        rdomains = []
        rchains = {}
        rfixed = []
        rsymtemp = None
        rmulti = []

        for j in range(len(args.chain[i])):
            # For each component of the chain

            if args.chain[i][j][0][-4:]=='.pdb':
                pdb = b.PDBModel(args.chain[i][j][0])
                    if len(args.chain[i][j])==2:  
                        # If chain to be taken is specified
                        rchains[pdb] = args.chain[i][j][1]
                        # CHANGE THIS TO WORK WITH DUPLICATED PDB NAMES
                        # MAYBE AN INDEX INSTEAD ?
                        # SYMTEMPLATE PDB CANNOT BE DUPLICATED
                        rmulti.append((args.chain[i][j][0], 
                            args.chain[i][j][1]))
                    if args.chain[i][j][0] in args.fixed:  
                        # If domain will be fixed
                        rfixed.append(pdb)
                    if args.chain[i][j][0] in args.symtemplate:  
                        # If is symtemplate ... THERE CAN ONLY BE ONE
                        rsymtemp = pdb

                rdomains.append(pdb)
                
            else:
                rdomains.append(args.chain[i][j][0])

        # chains[i] = (domains/linkers, chains dict, fixed domains, 
        #               symtemplate, pool symmetry, already modeled, 
        #               multichain domains)
        # rmulti is a list of tuples [('1234.pdb', 'A'), ... ] for each domain with
        # a chain specified
        chains.append((rdomains, rchains, rfixed, rsymtemp, args.poolsym, 
            False, rmulti))

    random.shuffle(chains)
    return chains

## IF NAME=__MAIN__

chains = create_args(args)

# Create models
# returns list with 10 models
call = builder.Builder(chains, args.destination, args.debug)
models = call.build()


#############
##  TESTING        
#############
import testing

class TestMultiprot(testing.AutoTest):
    """
    Test class

    Test for the same examples as in ranch.py (1, 4, 5, 7, 10)
    ONLY SINGLE CHAINS FOR NOW
    """

    def setUp(self):
        linker = linker or 'G'*15
        self.argstring1 = '--chain dom1.pdb '+linker+' dom2.pdb'
        self.argstring4 = '--chain domAB1.pdb:A '+linker+' domAB1.pdb:B'
        self.argstring5 = '--chain domAB1.pdb '+linker+' domAB2:A \
        --symmetry p2 --symtemplate domAB1.pdb --poolsym s'
        self.argstring7 = '--chain histone.pdb '+linker+' domAB1.pdb '+linker+
            ' histone.pdb --symmetry p2 --symtemplate domAB1.pdb --poolsym mix'
        self.argstring10 = '--chain domAB1.pdb '+linker+' domAB1.pdb:B '\
            +linker+' domAB1.pdb:B'

    def test_parsing(self):


