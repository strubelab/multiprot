#############
##  TESTING        
#############
import os
import biskit as B
import multiprot_script as mp
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
        args = mp.parsing(self.argstring1.split())
        CHAINS = mp.create_chains(args)

        self.assertTrue(len(CHAINS)==1, 'Incorrect number of chains')
        
        chain = CHAINS[0]
        self.assertTrue(isinstance(chain,mp.Chain))
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
        args = mp.parsing(self.argstring4.split())
        CHAINS = mp.create_chains(args)

        self.assertTrue(len(CHAINS)==1, 'Incorrect number of chains')
        
        chain = CHAINS[0]
        self.assertTrue(isinstance(chain,mp.Chain))
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
        args = mp.parsing(self.argstring5.split())
        CHAINS = mp.create_chains(args)

        self.assertTrue(len(CHAINS)==1, 'Incorrect number of chains')
        
        chain = CHAINS[0]

        self.assertTrue(isinstance(chain,mp.Chain))
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
        args = mp.parsing(self.argstring2ch.split())
        CHAINS = mp.create_chains(args)

        self.assertTrue(len(CHAINS)==2, 'Incorrect number of chains')

        chain0 = CHAINS[0]
        self.assertTrue(isinstance(chain0,mp.Chain))
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
        self.assertTrue(isinstance(chain1,mp.Chain))
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
        args = mp.parsing(self.argstring2chfixed.split())
        CHAINS = mp.create_chains(args)

        self.assertTrue(len(CHAINS)==2, 'Incorrect number of chains')
        
        chain0 = CHAINS[0]
        self.assertTrue(isinstance(chain0,mp.Chain))
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
        self.assertTrue(isinstance(chain1,mp.Chain))
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
        args = mp.parsing(self.argstring3ch.split())
        CHAINS = mp.create_chains(args)

        self.assertTrue(len(CHAINS)==3, 'Incorrect number of chains')
        
        chain0 = CHAINS[0]
        self.assertTrue(isinstance(chain0,mp.Chain))
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
        self.assertTrue(isinstance(chain1,mp.Chain))
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
        self.assertTrue(isinstance(chain2,mp.Chain))
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
