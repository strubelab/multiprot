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


