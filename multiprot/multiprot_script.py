"""
Script that handles command line argument parsing and calls Builder module

"""

import argparse
import biskit as b
import random
import os

import multiprot.builder as builder

# If type=divide does not work, try action='append_const'

def divide(s):
    """
    Converts entries for the --chain argument into tuples, and checks if PDB file
    exists

    E.g.
    ABCD.pdb:A --> ('ABCD.pdb','A')
    ABCD.pdb   --> ('ABCD.pdb',)
    
    :param s: string with one of the entries for --chain argument
    :type s: str

    """
    r = tuple(s.split(':'))
    # Does not work... find out why
    # if not os.path.exists(r[0]):
    #     raise argparse.ArgumentTypeError('Specified PDB file does not exist')
    
    return r

def path_exists(s):
    """
    Checks if path specified for output exists

    :param s: path specified for the --destination parameter
    :type s: str
    """
    if not os.path.exists(s):
        raise argparse.ArgumentTypeError('Invalid destination for output files.')
    return s


def parsing(args=None):
    """
    Creates the argument parser instance and applies it to the command line input

    :param args:    list with the arguments to be parsed (only for testing
                    purposes). If none is provided, it takes them from sys.argv
    :type args:     list
    """

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

    # Feature not supported yet
    # parser.add_argument('--split', '-spl', default=None, 
    #     help='Split one or more chains from a given PDB')

    parser.add_argument('--symmetry', '-sym', default='p1',
        help='What kind of symmetry do you wish to have in your molecule. Supported\n\
        symmetries are: p1, p2, ..., p19 (nineteen-fold), p22, p32, p42, p52, p62,\n\
        ..., p122, p222.', choices=supp_sym)

    parser.add_argument('--symtemplate', '-t', default=[], action='append', 
        help='Which domain will be the symmetry core, in case of symmetry other than\
         p1 specified')

    parser.add_argument('--poolsym', '-o', default='m', choices=['m', 's', 'a'], 
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

    return parser.parse_args(args)
    # vars(args)  returns dictionary with attributes


def create_chains(args):
    """
    Function that takes the arguments parsed from the command line and
    returns a list with all input necessary for ranch for each chain

    :param args:    Namespace object created by calling .parse_args() method on 
                    the given arguments (generally from sys.argv) 
    """
    
    chains = []

    for chain in args.chain:    # For each chain
        rdomains = []   # List of PDBModels and strings composing the chain
        rchains = {}    # Dictionary with chain specification for multichain pdbs
        rfixed = []     # List with domains to be fixed in their coordinates
        rsymtemp = None     # symtemplate
        rmulti = []     # list of tuples [('i_1234.pdb', 'A'), ... ] for each domain
                        # with a chain specified, where i is the domain number
                        # inside the chain

        for i in range(len(chain)):
            # For each component of the chain

            if chain[i][0][-4:]=='.pdb':     # If the element is a pdb structure
                pdb = b.PDBModel(chain[i][0])
                if len(chain[i])==2:  
                    # If the chain to be taken from the domain is specified, e.g.
                    # # ABCD.pdb:A --> ('ABCD.pdb','A')
                    rchains[pdb] = chain[i][1]
                    
                    # Append domain index to deal with repeated pdb names
                    rmulti.append((str(i)+'_'+chain[i][0], chain[i][1]))
                    # CHANGE THIS TO WORK WITH DUPLICATED PDB NAMES
                    # MAYBE APPEND AN INDEX ?
                    # SYMTEMPLATE PDB CANNOT BE DUPLICATED
                    
                if chain[i][0] in args.fixed:  
                    # If domain will be fixed
                    rfixed.append(pdb)
                
                if chain[i][0] in args.symtemplate:  
                    # If the domain is symtemplate ... THERE CAN ONLY BE ONE
                    rsymtemp = pdb

                rdomains.append(pdb)
                
            else:   # If the element is a linker (string with sequence of AA)
                rdomains.append(chain[i][0])

        # dictionary with arguments to be passed on to ranch
        args_dict = {
            "chains" : rchains, 
            "symmetry" : args.symmetry,
            "symtemplate" : rsymtemp, 
            "pool_sym" : args.poolsym,
            "fixed" : rfixed 
            }
        
        chains.append((rdomains, args_dict, False, rmulti))
        # chains[i] = (domains/linkers, args_dict, already modeled (T/F), 
        #               multichain domains)
        # rmulti is a list of tuples [('1234.pdb', 'A'), ... ] for every domain 
        # with a chain specified
        
    random.shuffle(chains)
    return chains


if __name__ == '__main__':

    # Parse arguments
    args = parsing()  
    chains = create_chains(args)

    # Create models
    # returns list with 10 models
    call = builder.Builder(chains, args.destination, args.debug)
    models = call.build()
