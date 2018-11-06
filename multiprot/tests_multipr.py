#############
##  TESTING        
#############
import tempfile, os
import multiprot.testing as testing
import multiprot.parseChains as C
import multiprot.builder as bu

class TestMultipr(testing.AutoTest):
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
    tetramer = None
    linker = None

    def setUp(self):

        self.testpath = self.testpath or \
            os.path.join(os.path.abspath(os.path.dirname(__file__)), 'testdata')

        self.mono1 = self.mono1 or os.path.join(self.testpath, '2z6o.pdb')
        self.mono2 = self.mono2 or os.path.join(self.testpath, 'histone.pdb')
        self.mono3 = self.mono3 or os.path.join(self.testpath, '1it2_A.pdb')
        self.dimer1 = self.dimer1 or os.path.join(self.testpath, 'domAB1.pdb')
        self.dimer2 = self.dimer2 or os.path.join(self.testpath, 'domAB2.pdb')
        self.dimer3 = self.dimer3 or os.path.join(self.testpath, '2qud_mod.pdb')
        self.trimer = self.trimer or os.path.join(self.testpath, '2ei4_mod.pdb')
        self.tetramer = self.tetramer or os.path.join(self.testpath, '5agc.pdb')

        self.linker = self.linker or 'TG'*25

        ## ADD TEST WITH 3 CHAINS AND SYMMETRY
    
    # PASSED
    def test_example1(self):
        '''
        Single chain example
        '''
        # testdir = tempfile.mkdtemp('', self.__class__.__name__.lower() + \
        #     '_example1_')

        argstring = '--chain '+self.mono1+' '+self.linker+' '+self.mono2# +\
            # ' --destination '+testdir ... to write the models to disk

        args = C.parsing(argstring.split())
        CHAINS = C.create_chains(args)
        build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

        model = build.run()

        self.assertTrue(model.lenChains()==1)
        #self.assertTrue(len(model)==2336)

        # build.write_pdbs([model],testdir)

    # PASSED
    def test_example4(self):
        '''
        Single chain with two multiple-chain domains
        '''
        
        argstring = \
            '--chain '+self.dimer1+':A '+self.linker+' '+self.dimer2+':B '+\
            '--fixed '+self.dimer1

        args = C.parsing(argstring.split())
        CHAINS = C.create_chains(args)
        build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

        model = build.run()

        self.assertTrue(model.lenChains()==3)
        #self.assertTrue(len(model)==7326)


    # # PASSED
    # def test_example5(self):
    #     '''
    #     Single chain with symmetry
    #     '''

    #     argstring = \
    #         '--chain '+self.dimer1+' '+self.linker+' '+self.dimer2+':A '+\
    #         '--symmetry p2 --symtemplate '+self.dimer1+' --poolsym s'

    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==4)
    #     #self.assertTrue(len(model)==11080)


    # FAILED
    # def test_2Chainz(self):
    #     '''
    #     TWO CHAAAAINZ
    #     This one might fail, because when the domains to be bound are not fixed 
    #     in their coordinates it is less likely that the program will be able to 
    #     connect them in the second chain, i.e. the model of the first chain might
    #     take them too far apart from each other
    #     '''

    #     argstring = \
    #         '--chain '+self.dimer1+':A '+self.linker+' '+self.dimer3+':A '+\
    #         self.linker+' '+self.dimer2+':A '+\
    #         '--chain '+self.dimer1+':B '+self.linker+' '+self.dimer3+':B '+\
    #         self.linker+' '+self.dimer2+':B'

    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==2)



    # # PASSED
    # def test_2Chainzfixed(self):
    #     '''
    #     Same as 2Chainz but the domains are fixed in their coordinates
    #     '''

    #     argstring = \
    #         '--chain '+self.dimer1+':A '+self.linker+' '+self.dimer3+':A '+\
    #         self.linker+' '+self.dimer2+':A '+\
    #         '--chain '+self.dimer1+':B '+self.linker+' '+self.dimer3+':B '+\
    #         self.linker+' '+self.dimer2+':B '+\
    #         '--fixed '+self.dimer1+' '+self.dimer2+' '+self.dimer3

    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==2)
    #     self.assertTrue(len(model)==10140)


    # PASSED
    # def test_symp2(self):
    #     '''
    #     Test with symmetry
    #     '''

    #     # There is no need to specify chain for the symmetric core
    #     argstring = \
    #         '--chain '+self.dimer1+' '+self.linker+' '+self.mono1+\
    #         ' --symmetry p2 --symtemplate '+self.dimer1+\
    #         ' --poolsym s'


    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==2)


    # PASSED
    def test_symp3(self):
        '''
        Test with symmetry
        '''

        # There is no need to specify chain for the symmetric core
        argstring = \
            '--chain '+self.trimer+' '+self.linker+' '+self.mono1+\
            ' --symmetry p3 --symtemplate '+self.trimer+\
            ' --poolsym s'


        args = C.parsing(argstring.split())
        CHAINS = C.create_chains(args)
        build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

        model = build.run()

        self.assertTrue(model.lenChains()==3)


    # NEED STRUCTURE WITH CHAINS EXACTLY EQUAL
    # FAILED
    # def test_symp4(self):
    #     '''
    #     Test with symmetry
    #     '''

    #     # There is no need to specify chain for the symmetric core
    #     argstring = \
    #         '--chain '+self.tetramer+' '+self.linker+' '+self.mono1+\
    #         ' --symmetry p4 --symtemplate '+self.tetramer+\
    #         ' --poolsym s'


    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==4)



    # PASSED
    # def test_4ch1(self):
    #     '''
    #     Example 1 with three chains, using a trimer
    #     '''

    #     argstring = \
    #         '--chain '+self.mono1+' '+self.linker+' '+self.tetramer+':A '+\
    #         self.linker+' '+self.mono1+\
    #         ' --chain '+self.mono2+' '+self.linker+' '+self.tetramer+':B '+\
    #         self.linker+' '+self.mono2+\
    #         ' --chain '+self.mono3+' '+self.linker+' '+self.tetramer+':C '+\
    #         self.linker+' '+self.mono3+\
    #         ' --chain '+self.mono2+' '+self.linker+' '+self.tetramer+':D '+\
    #         self.linker+' '+self.mono3

    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==4)


    # PASSED
    # def test_4ch2(self):
    #     '''
    #     Example 1 with three chains, using a trimer
    #     '''

    #     argstring = \
    #         '--chain '+self.tetramer+':A '+\
    #         self.linker+' '+self.mono1+\
    #         ' --chain '+self.tetramer+':B '+\
    #         self.linker+' '+self.mono2+\
    #         ' --chain '+self.tetramer+':C '+\
    #         self.linker+' '+self.mono3+\
    #         ' --chain '+self.tetramer+':D '+\
    #         self.linker+' '+self.mono3

    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==4)


    # def test_2Chainzsym(self):
    #     '''
    #     Test with two chains and three-fold symmetry
    #     '''
    #     tempdir = tempfile.mkdtemp('', self.__class__.__name__.lower() + \
    #         '_2chsym_')

    #     # There is no need to specify chain for the symmetric core
    #     argstring = \
    #         '--chain '+self.trimer+' '+self.linker+' '+self.dimer1+':A '+\
    #         '--chain '+self.mono1+' '+self.linker+' '+self.dimer1+':B '+\
    #         self.linker+' '+self.mono2+\
    #         ' --symmetry p3 --symtemplate '+self.trimer+\
    #         ' --poolsym s'


    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()
    #     # build.write_pdbs([model],tempdir)
    #     self.assertTrue(model.lenChains()==6)


    # # PASSED
    # def test_3ch1(self):
    #     '''
    #     Example 1 with three chains, using a trimer
    #     '''

    #     argstring = \
    #         '--chain '+self.trimer+':A '+self.linker+' '+self.mono1+\
    #         ' --chain '+self.trimer+':B '+self.linker+' '+self.mono2+\
    #         ' --chain '+self.trimer+':C '+self.linker+' '+self.mono3

    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==3)


    # PASSED
    def test_3ch2(self):
        '''
        Example 2 with three chains
        '''

        argstring = \
            '--chain '+self.mono1+' '+self.linker+' '+self.dimer1+':A '+\
            self.linker+' '+self.mono2+\
            ' --chain '+self.dimer1+':B '+self.linker+' '+self.dimer2+':A'+\
            ' --chain '+self.mono1+' '+self.linker+' '+self.dimer2+':B '+\
            self.linker+' '+self.mono2

        args = C.parsing(argstring.split())
        CHAINS = C.create_chains(args)
        build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

        model = build.run()

        self.assertTrue(model.lenChains()==3)


    # # PASSED
    # def test_3ch3(self):
    #     '''
    #     Example 3 with three chains, using a trimer
    #     '''

    #     argstring = \
    #         '--chain '+self.mono1+' '+self.linker+' '+self.trimer+':A '+\
    #         self.linker+' '+self.mono1+\
    #         ' --chain '+self.mono1+' '+self.linker+' '+self.trimer+':B '+\
    #         self.linker+' '+self.mono2+\
    #         ' --chain '+self.mono1+' '+self.linker+' '+self.trimer+':C '+\
    #         self.linker+' '+self.mono3+\
    #         ' --destination '+tempdir

    #     args = C.parsing(argstring.split())
    #     CHAINS = C.create_chains(args)
    #     build = bu.Builder(CHAINS,args.debug,args.number,args.destination)

    #     model = build.run()

    #     self.assertTrue(model.lenChains()==3)


if __name__ == '__main__':

    testing.localTest(debug=False)

