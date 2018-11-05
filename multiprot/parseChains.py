#!/usr/bin/env python3
"""
Script that handles command line argument parsing and calls Builder module

"""

import argparse
import biskit as B
import random
import os

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
        os.makedirs(s)
    return s

def number_models(n):
    if n<10:
        return 10
    else:
        return n

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

    parser = argparse.ArgumentParser(prog='multipr', 
        description='''Build protein models connecting two or more structured 
        domains with disordered linkers. Supports modelling of multiple chains 
        and symmetry.''')

    parser.add_argument('--chain', '-c', action='append', nargs='+', type=divide, 
        help='Add a new chain to the model', required=True)

    # Feature not supported yet
    # parser.add_argument('--split', '-spl', default=None, 
    #     help='Split one or more chains from a given PDB')

    parser.add_argument('--symmetry', '-sym', default='p1',
        help='What kind of symmetry do you wish to have in your molecule. Supported\n\
        symmetries are: p1, p2, ..., p19 (nineteen-fold), p22, p32, p42, p52, p62,\n\
        ..., p122, p222.')  # choices=supp_sym

    parser.add_argument('--symtemplate', '-t', default=[], 
        help='Which domain will be the symmetry core, in case of symmetry other than\
         p1 specified')

    parser.add_argument('--number', '-n', default=3, help='How many models do you want\
        to produce? (less models = faster)')

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
    Function that takes the arguments parsed from the command line and returns
    a list of Chain objects with all input necessary for ranch for each chain

    :param args:    Namespace object created by calling .parse_args() method on 
                    the given arguments (generally from sys.argv) 
    """
    
    CHAINS = []     # List of Chain objects, one for each chain in the model

    for chain in args.chain:    # For each chain
        rnames = []     # Arguments for --chain as provided in input
        rdomains = []   # List of PDBModels and strings composing the chain
        rchains = {}    # Dictionary with chain specification for multichain pdbs
        rchains_names = {}  # Same as rchains above, but the keys are names instead
                            # of PDBModels
        rfixed = []     # List with domains to be fixed in their coordinates
        rsymtemp = None     # symtemplate
        # rmulti = []     # list of tuples [('i_1234.pdb', 'A'), ... ] for each domain
                        # with a chain specified, where i is the domain number
                        # inside the chain

        for i in range(len(chain)):
            # For each component of the chain
            rnames.append(chain[i])
            if chain[i][0][-4:]=='.pdb':     # If the element is a pdb structure
                pdb = B.PDBModel(chain[i][0])
                if len(chain[i])==2:  
                    # If the chain to be taken from the domain is specified, e.g.
                    # # ABCD.pdb:A --> ('ABCD.pdb','A')
                    rchains[pdb] = chain[i][1]
                    rchains_names[chain[i][0]] = chain[i][1]
                    
                    # Append domain index to deal with repeated pdb names
                    # rmulti.append((str(i)+'_'+chain[i][0], chain[i][1]))
                    # CHANGE THIS TO WORK WITH DUPLICATED PDB NAMES
                    # MAYBE APPEND AN INDEX ?
                    # SYMTEMPLATE PDB CANNOT BE DUPLICATED
                    
                if chain[i][0] in args.fixed:  
                    # If domain will be fixed
                    rfixed.append(pdb)

                    # Remove the pdb name from args.fixed so it won't be duplicated
                    # if it is in multiple chains
                    # Should I keep a copy of the original args.fixed?
                    args.fixed.remove(chain[i][0])
                
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
            "fixed" : rfixed,
            "symunit" : None,
            "n" : number_models(args.number)
            }

        CHAINS.append(Chain(rnames, rdomains, args_dict, False, rchains_names))

        # chains[i] = (domains/linkers, args_dict, already modeled (T/F), 
        #               multichain domains)
        # rmulti is a list of tuples [('1234.pdb', 'A'), ... ] for every domain 
        # with a chain specified
        
    # Put any chain with fixed domains at the beginning, to be modeled first
    if any([ch.args["fixed"] for ch in CHAINS]):
        while not CHAINS[0].args["fixed"]:
            random.shuffle(CHAINS)
    # else:
        # For testing, it is necessary to remove shuffling of chains
        # random.shuffle(CHAINS)  # Necessary?

    return CHAINS


class Chain:
    """
    Do your chain hang low
    do it wobble to the floor
    do it shine in the light
    is it platinum is it gold
    """

    def __init__(self, names, domains, args, modeled, chains_names):
        # Super weird behaviour when putting these variabels in the init arguments
        # paired_to=[],modeled_domains=[],emb_seq=None,container_seq=None,
        # jdomains={},testvar=[]):
        self.names = names
        self.domains = domains
        self.args = args
        self.modeled = modeled
        self.chains_names = chains_names
        self.paired_to = []
        # modeled_domains of every symmetric unit
        self.modeled_domains = []
        self.emb_seq = None
        self.container_seq = None

        # Domains with new coordinates for non-symmetric structures
        self.new_domains = ['']*len(self.names)
        # self.testvar= testvar
        
        # Dictionary with the modeled symmetric domains of the current chain. 
        # Similar to modeled_domains, but this only contains the
        # part of the domain belonging to this chain, i.e. what was removed from
        # self.full_chains at the end of builder.replace_modeled()
        # These domains will be used to embed the symmetric units in builder.run()
        # The dictionary will be of the form {name_index:[domains]}, where
        # nameindex is the place or order of the domain in the chain, and [domains]
        # is a list with the symmetric domains' coordinates (only 1 if there is no
        # symmetry)
        self.jdomains = {}


#############
##  TESTING        
#############
import os
import biskit as B
import multiprot.testing as testing

class TestMultiprot(testing.AutoTest):
    """
    Test class for argument parsing

    Test for the examples 1, 4, 5 and 
    ONLY SINGLE CHAINS FOR NOW
    """

    testpath = None
    mono1 = None
    mono2 = None
    mono3 = None
    dimer1 = None
    dimer2 = None
    dimer3 = None
    trimer = None
    linker = None
    argstring1 = None
    argstring4 = None
    argstring5 = None
    argstring2ch = None
    argstring2chfixed = None
    argstring3ch = None

    def setUp(self):

        self.testpath = self.testpath or \
            os.path.join(os.path.abspath(os.path.dirname(__file__)), 'testdata')

        self.mono1 = self.mono1 or os.path.join(self.testpath, '2z6o.pdb')
        self.mono2 = self.mono2 or os.path.join(self.testpath, 'histone.pdb')
        self.mono3 = self.mono3 or os.path.join(self.testpath, '2h5q.pdb')
        self.dimer1 = self.dimer1 or os.path.join(self.testpath, 'domAB1.pdb')
        self.dimer2 = self.dimer2 or os.path.join(self.testpath, 'domAB2.pdb')
        self.dimer3 = self.dimer3 or os.path.join(self.testpath, '2qud.pdb')
        self.trimer = self.trimer or os.path.join(self.testpath, '2ei4.pdb')

        self.linker = self.linker or 'TG'*15

        #### SINGLE CHAIN EXAMPLES
        # Assemble argument strings
        self.argstring1 = self.argstring1 or \
                '--chain '+self.mono1+' '+self.linker+' '+self.mono2
        
        # Careful with the names of fixed domains... cannot be repeated in the chain
        # for repeating pdbs that are not fixed or are not symtemplates there
        # should be no problem
        self.argstring4 = self.argstring4 or \
                '--chain '+self.dimer1+':A '+self.linker+' '+\
                self.dimer2+':B --fixed '+self.dimer1
        
        #### SYMMETRY EXAMPLE
        # Careful with symtemplate name... cannot be repeated in the chain
        self.argstring5 = self.argstring5 or '--chain '+self.dimer1+' '+\
                self.linker+' '+self.dimer2+':A --symmetry p2 --symtemplate '\
                +self.dimer1+' --poolsym s'

        #### DOUBLE CHAIN
        # Not fixed
        self.argstring2ch = self.argstring2ch or '--chain '+self.dimer1+':A '+\
                self.linker+' '+self.dimer3+':A '+self.linker+' '+self.dimer2+':A '+\
                '--chain '+self.dimer1+':B '+self.linker+' '+self.dimer3+':B '+\
                self.linker+' '+self.dimer2+':B'

        # Fixed
        self.argstring2chfixed = self.argstring2chfixed or '--chain '+self.dimer1+\
                ':A '+self.linker+' '+self.dimer3+':A '+self.linker+' '+\
                self.dimer2+':A '+'--chain '+self.dimer1+':B '+self.linker+' '+\
                self.dimer3+':B '+self.linker+' '+self.dimer2+':B --fixed '+\
                self.dimer1+' '+self.dimer2+' '+self.dimer3

        #### TRIPLE CHAIN
        self.argstring3ch = self.argstring3ch or '--chain '+self.trimer+':A '+\
                self.linker+' '+self.mono1+' --chain '+self.trimer+':B '+\
                self.linker+' '+self.mono2+' --chain '+self.trimer+':C '+\
                self.linker+' '+self.mono3


    def test_example1(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'CHAINS' result from create_chains()
        """
        args = parsing(self.argstring1.split())
        CHAINS = create_chains(args)

        self.assertTrue(len(CHAINS)==1, 'Incorrect number of chains')
        
        chain = CHAINS[0]
        self.assertTrue(isinstance(chain,Chain))
        self.assertTrue(len(chain.names)==3)
        domains = chain.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel), 'Problem with domains list')
    
    def test_example4(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'CHAINS' result from create_chains()
        """
        args = parsing(self.argstring4.split())
        CHAINS = create_chains(args)

        self.assertTrue(len(CHAINS)==1, 'Incorrect number of chains')
        
        chain = CHAINS[0]
        self.assertTrue(isinstance(chain,Chain))
        self.assertTrue(len(chain.names)==3)
        domains = chain.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain.args["chains"].items()]))
        ch_val = [v for k,v in chain.args["chains"].items()]
        self.assertTrue(len(ch_val)==2 and ('A' in ch_val) and ('B' in ch_val),
            'Problem with chains dictionary')
        chnames = [(k,v) for k,v in chain.chains_names.items()]
        self.assertTrue(chnames == [(self.dimer1,'A'),(self.dimer2,'B')])
        self.assertTrue(len(chain.args["fixed"])==1 and \
            isinstance(chain.args["fixed"][0], B.PDBModel), "Problem with 'fixed' \
            argument")

    def test_example5(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'CHAINS' result from create_chains()
        """
        args = parsing(self.argstring5.split())
        CHAINS = create_chains(args)

        self.assertTrue(len(CHAINS)==1, 'Incorrect number of chains')
        
        chain = CHAINS[0]

        self.assertTrue(isinstance(chain,Chain))
        self.assertTrue(len(chain.names)==3)
        domains = chain.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain.args["chains"].items()]))
        ch_val = [v for k,v in chain.args["chains"].items()]
        self.assertTrue(ch_val == ['A'], 'Problem with chains dictionary')
        chnames = [(k,v) for k,v in chain.chains_names.items()]
        self.assertTrue(chnames == [(self.dimer2,'A')])
        self.assertTrue(chain.args["symmetry"]=="p2")
        self.assertTrue(isinstance(chain.args["symtemplate"], B.PDBModel))
        self.assertTrue(chain.args["pool_sym"]=="s")
        
    def test_2Chainz(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'CHAINS' result from create_chains()
        """
        args = parsing(self.argstring2ch.split())
        CHAINS = create_chains(args)

        self.assertTrue(len(CHAINS)==2, 'Incorrect number of chains')

        chain0 = CHAINS[0]
        self.assertTrue(isinstance(chain0,Chain))
        self.assertTrue(len(chain0.names)==5)
        domains = chain0.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel) and
            isinstance(domains[3], str) and 
            isinstance(domains[4], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain0.args["chains"].items()]))
        ch_val = [v for k,v in chain0.args["chains"].items()]
        self.assertTrue(ch_val == ['A', 'A', 'A'], 'Problem with chains dictionary')
        chnames = [(k,v) for k,v in chain0.chains_names.items()]
        self.assertTrue(chnames == [(self.dimer1,'A'),(self.dimer3,'A'),(self.dimer2,'A')])

        chain1 = CHAINS[1]
        self.assertTrue(isinstance(chain1,Chain))
        self.assertTrue(len(chain1.names)==5)
        domains = chain1.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel) and
            isinstance(domains[3], str) and 
            isinstance(domains[4], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain1.args["chains"].items()]))
        ch_val = [v for k,v in chain1.args["chains"].items()]
        self.assertTrue(ch_val == ['B', 'B', 'B'], 'Problem with chains dictionary')
        chnames = [(k,v) for k,v in chain1.chains_names.items()]
        self.assertTrue(chnames == [(self.dimer1,'B'),(self.dimer3,'B'),(self.dimer2,'B')])

    def test_2Chainz_fixed(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'CHAINS' result from create_chains()
        """
        args = parsing(self.argstring2chfixed.split())
        CHAINS = create_chains(args)

        self.assertTrue(len(CHAINS)==2, 'Incorrect number of chains')
        
        chain0 = CHAINS[0]
        self.assertTrue(isinstance(chain0,Chain))
        self.assertTrue(len(chain0.names)==5)
        domains = chain0.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel) and
            isinstance(domains[3], str) and 
            isinstance(domains[4], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain0.args["chains"].items()]))
        ch_val = [v for k,v in chain0.args["chains"].items()]
        self.assertTrue(ch_val == ['A', 'A', 'A'], 'Problem with chains dictionary')
        self.assertTrue(len(chain0.args["fixed"])==3 and \
            all([isinstance(p, B.PDBModel) for p in chain0.args["fixed"]]),
                "Problem with 'fixed' argument")
        chnames = [(k,v) for k,v in chain0.chains_names.items()]
        self.assertTrue(chnames == [(self.dimer1,'A'),(self.dimer3,'A'),(self.dimer2,'A')])

        chain1 = CHAINS[1]
        self.assertTrue(isinstance(chain1,Chain))
        self.assertTrue(len(chain1.names)==5)
        domains = chain1.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel) and
            isinstance(domains[3], str) and 
            isinstance(domains[4], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain1.args["chains"].items()]))
        ch_val = [v for k,v in chain1.args["chains"].items()]
        self.assertTrue(ch_val == ['B', 'B', 'B'], 'Problem with chains dictionary')
        self.assertTrue(len(chain1.args["fixed"])==0, "Problem with 'fixed' \
            argument")
        chnames = [(k,v) for k,v in chain1.chains_names.items()]
        self.assertTrue(chnames == [(self.dimer1,'B'),(self.dimer3,'B'),(self.dimer2,'B')])

    def test_3Chains(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'CHAINS' result from create_chains()
        """
        args = parsing(self.argstring3ch.split())
        CHAINS = create_chains(args)

        self.assertTrue(len(CHAINS)==3, 'Incorrect number of chains')
        
        chain0 = CHAINS[0]
        self.assertTrue(isinstance(chain0,Chain))
        self.assertTrue(len(chain0.names)==3)
        domains = chain0.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain0.args["chains"].items()]))
        ch_val = [v for k,v in chain0.args["chains"].items()]
        self.assertTrue(ch_val == ['A'], 'Problem with chains dictionary')
        chnames = [(k,v) for k,v in chain0.chains_names.items()]
        self.assertTrue(chnames == [(self.trimer,'A')])

        chain1 = CHAINS[1]
        self.assertTrue(isinstance(chain1,Chain))
        self.assertTrue(len(chain1.names)==3)
        domains = chain1.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain1.args["chains"].items()]))
        ch_val = [v for k,v in chain1.args["chains"].items()]
        self.assertTrue(ch_val == ['B'], 'Problem with chains dictionary')
        chnames = [(k,v) for k,v in chain1.chains_names.items()]
        self.assertTrue(chnames == [(self.trimer,'B')])

        chain2 = CHAINS[2]
        self.assertTrue(isinstance(chain2,Chain))
        self.assertTrue(len(chain2.names)==3)
        domains = chain2.domains
        self.assertTrue(isinstance(domains[0], B.PDBModel) and 
            isinstance(domains[1], str) and 
            isinstance(domains[2], B.PDBModel), 'Problem with domains list')
        self.assertTrue(all([isinstance(k,B.PDBModel) for k,v in \
            chain2.args["chains"].items()]))
        ch_val = [v for k,v in chain2.args["chains"].items()]
        self.assertTrue(ch_val == ['C'], 'Problem with chains dictionary')
        chnames = [(k,v) for k,v in chain2.chains_names.items()]
        self.assertTrue(chnames == [(self.trimer,'C')])


if __name__ == '__main__':

    testing.localTest(debug=False)

