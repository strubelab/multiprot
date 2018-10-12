"""
Script that contains the higher level implementation of the ranch wrapper,
to handle single and multiple chain scenarios

"""


##### FOR THE EARLY IMPLEMENTATION OF MULTIPROT, WE WILL ASSUME ONLY ONE CHAIN
##### IN THE INPUT

import tempfile
import os, re
import biskit as B
import biskit.tools as T
import numpy as N
import multiprot.ranch as R
import multiprot.pulchra as P

class Builder:
    """
    Class that will contain all the necessary methods to call the Ranch and
    pulchra wrappers multiple times and build single and multiple chain 
    structures. This class will be used by the multiprot script.
    
    """

    def __init__(self, chains, dest, debug):
        """
        :param args: Object that contains the arguments parsed from the command line
        :type args: argparse.Namespace object created by calling parser.parse_args()
        """
        self.chains = chains
        self.dest = dest
        self.debug = debug

    def cleanup(self):
        """
        Delete temporary files
        """
        if not self.debug:
            T.tryRemove(self.tempdir, tree=True)

    def build(self):
        """
        Method that will build the models and return a list with 10
        PDBModels
        """
        try:
            # Call ranch for single chain ... for there is only one chain
            call = R.Ranch(*self.chains[0][0], **self.chains[0][1], debug=self.debug)
            models = call.run()

            # create temporary folder to write models
            self.tempdir = tempfile.mkdtemp('', self.__class__.__name__.lower() + '_')
            # create names to write pdbs into
            pdb_paths = [os.path.join(self.tempdir, 'm%02d' % i) for i in range(1,11)]

            # Determine number of symunits
            n = int(self.chains[0][1]['symmetry'][1:])
            
            for i in range(len(models)):

                # Rebuild only the chains with CA, which is the first one plus its
                # repetitions in symmetric units

                full = B.PDBModel()
                m = models[i]
                rb_seq = m.takeChains([0]).sequence()
                rb_len = len(m.takeChains([0]))
                chains_not_rb = list(range(m.lenChains()))

                for j in range(m.lenChains()):
                    chain = m.takeChains([j])

                    if chain.sequence() == rb_seq and len(chain) == rb_len:
                        chname = pdb_paths[i]+'_ch'+str(j)+'.pdb'
                        chain.writePdb(chname)

                        call = P.Pulchra(chname)
                        rebuild = call.run()
                        full = full.concat(B.PDBModel(chname[:-3]+'rebuilt.pdb'))

                    else:
                        full = full.concat(chain)
            
                full.addChainId()
                full['serial_number'] = N.arange(1,len(full)+1)
                full.writePdb(os.path.join(self.dest, 'm%02d.pdb' % (i+1)))
        except:
            raise
        finally:
            self.cleanup()


#############
##  TESTING        
#############
import multiprot.testing as testing
import multiprot.multiprot_script as mp

class TestBuilder(testing.AutoTest):
    """
    Test class for model building

    Test for the examples 1, 4 and 5
    ONLY SINGLE CHAINS FOR NOW
    """

    testpath = None
    testdir = None
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

        self.testdir = tempfile.mkdtemp('', self.__class__.__name__.lower() + '_')

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

    def tearDown(self):
        T.tryRemove(self.testdir, tree=True)

    def test_example1(self):
        """ 
        Test model building with self.argstring1 and check the chains of the
        resulting model
        """
        args = mp.parsing(self.argstring1.split())
        chains = mp.create_chains(args)

        call = Builder(chains, self.testdir, False)
        build = call.build()

        f_out = [os.path.join(self.testdir, 'm%02d.pdb' % i) for i in range(1,11)]
        self.assertTrue(all([os.path.exists(f) for f in f_out]), 'Models not written.')

        m = B.PDBModel(f_out[0])
        self.assertTrue(len(m)==2231, 'Incorrect model length')
        self.assertTrue(m.lenChains()==1)

    def test_example4(self):
        args = mp.parsing(self.argstring4.split())
        chains = mp.create_chains(args)

        call = Builder(chains, self.testdir, False)
        build = call.build()

        f_out = [os.path.join(self.testdir, 'm%02d.pdb' % i) for i in range(1,11)]
        self.assertTrue(all([os.path.exists(f) for f in f_out]), 'Models not written.')

        m = B.PDBModel(f_out[0])
        self.assertTrue(len(m)==7221, 'Incorrect model length')
        self.assertTrue(m.lenChains()==3)

    def test_example5(self):
        args = mp.parsing(self.argstring5.split())
        chains = mp.create_chains(args)

        call = Builder(chains, self.testdir, False)
        build = call.build()

        f_out = [os.path.join(self.testdir, 'm%02d.pdb' % i) for i in range(1,11)]
        self.assertTrue(all([os.path.exists(f) for f in f_out]), 'Models not written.')

        m = B.PDBModel(f_out[0])
        self.assertTrue(len(m)==10870, 'Incorrect model length')
        self.assertTrue(m.lenChains()==4)


if __name__ == '__main__':

    testing.localTest(debug=False)
