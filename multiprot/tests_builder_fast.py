#############
##  TESTING        
#############
import multiprot.parseChains as C
import multiprot.testing as testing
import biskit as B
import multiprot.builder as bu
import os, tempfile

class TestBuilder(testing.AutoTest):
    """
    Test class for argument parsing

    Test for the examples 1, 4, 5 and 
    ONLY SINGLE CHAINS FOR NOW
    """

    tempdir = None
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
    builder1 = None
    builder4 = None
    builder5 = None
    builder2ch = None
    builder2chfixed = None
    builder3ch = None

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
        # The model is the same as the previous one, but the coordinates of the
        # dimers will be fixed
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

        # Create Chain objects and Builder instance
        args1 = C.parsing(self.argstring1.split())
        CHAINS1 = C.create_chains(args1)
        self.builder1 = bu.Builder(CHAINS1,args1.debug,args1.number,
            args1.destination)

        args4 = C.parsing(self.argstring4.split())
        CHAINS4 = C.create_chains(args4)
        self.builder4 = bu.Builder(CHAINS4,args4.debug,args4.number,
            args4.destination)

        args5 = C.parsing(self.argstring5.split())
        CHAINS5 = C.create_chains(args5)
        self.builder5 = bu.Builder(CHAINS5,args5.debug,args5.number,
            args5.destination)

        args2ch = C.parsing(self.argstring2ch.split())
        CHAINS2ch = C.create_chains(args2ch)
        self.builder2ch = bu.Builder(CHAINS2ch,args2ch.debug,args2ch.number,
            args2ch.destination)

        args2chfixed = C.parsing(self.argstring2chfixed.split())
        CHAINS2chfixed = C.create_chains(args2chfixed)
        self.builder2chfixed = bu.Builder(CHAINS2chfixed,args2chfixed.debug,
            args2chfixed.number,args2chfixed.destination)

        args3ch = C.parsing(self.argstring3ch.split())
        CHAINS3ch = C.create_chains(args3ch)
        self.builder3ch = bu.Builder(CHAINS3ch,args3ch.debug,args3ch.number,
            args3ch.destination)

        ## ADD TEST WITH 3 CHAINS AND SYMMETRY
    
    # PASSED
    def test_find_paired(self):
        '''
        Test the output of find_paired method, which should return a dictionary 
        of the form 

        paired_to_i = {j:[pair_ij1,pair_ij2],k:[pair_ik1],...}

        Where 
        - i is the index of the chain whose bound chains will be found
        - j and k are indexes of bound chains
        - pair_ij1 and pair_ij2 are the names of the domains that bind chains i and j,
          in the form:
            pair_ijx = [(pdb_namei, chain_idi),(pdb_namej, chain_idj)]
        - pair_ik1 has the names of the domain that binds chains i and k, in the
          form:
            pair_ikx = [(pdb_namei, chain_idi),(pdb_namek, chain_idk)]

        This test method tests all the example cases in the setUp method
        
        '''

        # Only one chain with linkers, no chains paired
        paired_to1 = self.builder1.find_paired(0)
        self.assertTrue(len(paired_to1)==0)

        paired_to4 = self.builder4.find_paired(0)
        self.assertTrue(len(paired_to4)==0)

        # Symmetric chains, no paired chains taken into account
        paired_to5 = self.builder5.find_paired(0)
        self.assertTrue(len(paired_to5)==0)

        # Two chains with linkers
        # Running method on chain 0
        paired_to2ch0 = self.builder2ch.find_paired(0)
        pair1 = [(self.dimer1,'A'),(self.dimer1,'B')]
        pair2 = [(self.dimer3,'A'),(self.dimer3,'B')]
        pair3 = [(self.dimer2,'A'),(self.dimer2,'B')]
        self.assertTrue(paired_to2ch0=={1:[pair1,pair2,pair3]}, paired_to2ch0)
        # Running method on chain 1
        paired_to2ch1 = self.builder2ch.find_paired(1)
        pair1 = [(self.dimer1,'B'),(self.dimer1,'A')]
        pair2 = [(self.dimer3,'B'),(self.dimer3,'A')]
        pair3 = [(self.dimer2,'B'),(self.dimer2,'A')]
        self.assertTrue(paired_to2ch1=={0:[pair1,pair2,pair3]}, paired_to2ch1)

        # Three chains with linkers
        # Running method on chain 0
        paired_to3ch0 = self.builder3ch.find_paired(0)
        pair1 = [(self.trimer,'A'),(self.trimer,'B')]
        pair2 = [(self.trimer,'A'),(self.trimer,'C')]
        self.assertTrue(paired_to3ch0=={1:[pair1],2:[pair2]})
        # Running on chain 1
        paired_to3ch1 = self.builder3ch.find_paired(1)
        pair1 = [(self.trimer,'B'),(self.trimer,'A')]
        pair2 = [(self.trimer,'B'),(self.trimer,'C')]
        self.assertTrue(paired_to3ch1=={0:[pair1],2:[pair2]})
        # Running on chain 2
        paired_to3ch2 = self.builder3ch.find_paired(2)
        pair1 = [(self.trimer,'C'),(self.trimer,'A')]
        pair2 = [(self.trimer,'C'),(self.trimer,'B')]
        self.assertTrue(paired_to3ch2=={0:[pair1],1:[pair2]})

    # PASSED
    def test_embed_symmetric(self):
        '''
        Tests builder.embed_symmetric()
        '''
        mod1 = B.PDBModel(os.path.join(self.testpath, 'chain01_2ch.pdb'))
        emb_mod = mod1.takeChains([1,2,3])
        j_dom = mod1.takeChains([0])
        full_emb = self.builder1.embed_symmetric([j_dom],[emb_mod])

        self.assertTrue(len(full_emb[0])==9267, str(len(full_emb)))

    # PASSED  
    def test_extract_embedded(self):
        """
        Tests the builder.extract_embedded method
        """
        
        mod1 = B.PDBModel(os.path.join(self.testpath, 'chain01_2ch.pdb'))
        emb_mod = mod1.takeChains([1,2,3])
        j_dom = mod1.takeChains([0])
        full_emb = self.builder1.embed_symmetric([j_dom],[emb_mod])
        container_seq = j_dom.sequence()[:2] + emb_mod.sequence() +\
            j_dom.sequence()[2:]

        full = full_emb[0]
        while full.lenChains()>1:
            full.mergeChains(0)

        chain01_2ch_reb = self.builder1.extract_embedded(full, emb_mod, 
            container_seq)

        # chain01_2ch_reb.writePdb('testdata/chain01_testrebuilt.pdb')

        self.assertTrue(chain01_2ch_reb.lenChains()==4)
        self.assertTrue(chain01_2ch_reb.sequence()==mod1.sequence())
        # The length is 9747 instead of 9750 because Biskit does not write OXT
        # And pulchra keeps it only for the rebuilt chain
        self.assertTrue(len(chain01_2ch_reb)==9747, len(chain01_2ch_reb))
        


if __name__ == '__main__':

    testing.localTest(debug=False)
