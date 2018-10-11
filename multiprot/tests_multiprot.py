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

    Test for the examples 1, 4 and 5
    ONLY SINGLE CHAINS FOR NOW
    """

    testpath = None
    dom1path = None
    dom2path = None
    domAB1path = None
    linker = None
    argstring1 = None
    argstring4 = None
    argstring5 = None

    def setUp(self):

        self.testpath = self.testpath or \
            os.path.join(os.path.abspath(os.path.dirname(__file__)), 'testdata')

        self.dom1path = os.path.join(self.testpath, '2z6o.pdb')
        self.dom2path = os.path.join(self.testpath, 'histone.pdb')
        self.domAB1path = os.path.join(self.testpath, 'domAB1.pdb')

        self.linker = self.linker or 'G'*15

        # Assemble argument strings
        self.argstring1 = self.argstring1 or \
                '--chain '+self.dom1path+' '+self.linker+' '+self.dom2path
        
        self.argstring4 = self.argstring4 or \
                '--chain '+self.domAB1path+':A '+self.linker+' '+\
                self.domAB1path+':B'
        
        self.argstring5 = self.argstring5 or '--chain '+self.domAB1path+' '+\
                self.linker+' '+self.domAB1path+':A --symmetry p2 --symtemplate '\
                +self.domAB1path+' --poolsym s'
        
    def test_example1(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'chains' result from create_chains()
        """
        args = mp.parsing(self.argstring1.split())
        chains = mp.create_chains(args)

        self.assertTrue(len(chains)==1, 'Incorrect number of chains')
        chain = chains[0]

        self.assertTrue(len(chain)==4, 'Chain does not have all the elements')
        self.assertTrue(isinstance(chain[0][0], B.PDBModel) and 
            isinstance(chain[0][1], str) and 
            isinstance(chain[0][2], B.PDBModel), 'Problem with domains list')
    
    def test_example4(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'chains' result from create_chains()
        """
        args = mp.parsing(self.argstring4.split())
        chains = mp.create_chains(args)

        self.assertTrue(len(chains)==1, 'Incorrect number of chains')
        chain = chains[0]

        self.assertTrue(len(chain)==4, 'Chain does not have all the elements')
        self.assertTrue(isinstance(chain[0][0], B.PDBModel) and 
            isinstance(chain[0][1], str) and 
            isinstance(chain[0][2], B.PDBModel), 'Problem with domains list')
        self.assertTrue(len(chain[1]["chains"])==2 and all([isinstance(v, str) 
            for k,v in chain[1]["chains"].items()]), 
            'Problem with chains dictionary')

    def test_example5(self):
        """ 
        Test argument parsing from self.argstring1 and check the elements of the 
        'chains' result from create_chains()
        """
        args = mp.parsing(self.argstring5.split())
        chains = mp.create_chains(args)

        self.assertTrue(len(chains)==1, 'Incorrect number of chains')
        chain = chains[0]

        self.assertTrue(len(chain)==4, 'Chain does not have all the elements')
        self.assertTrue(isinstance(chain[0][0], B.PDBModel) and 
            isinstance(chain[0][1], str) and 
            isinstance(chain[0][2], B.PDBModel), 'Problem with domains list')
        self.assertTrue(len(chain[1]["chains"])==1 and all([isinstance(v, str) 
            for k,v in chain[1]["chains"].items()]), 
            'Problem with chains dictionary')
        self.assertTrue(chain[1]["symmetry"]=="p2")
        self.assertTrue(isinstance(chain[1]["symtemplate"], B.PDBModel))
        self.assertTrue(chain[1]["pool_sym"]=="s")

if __name__ == '__main__':

    testing.localTest(debug=False)
