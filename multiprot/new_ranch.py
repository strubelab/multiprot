"""
Author: Francisco Javier Guzman

Last edited: 20/feb/18

"""


import biskit as B
import numpy as N
import tempfile, os

from biskit.exe.executor import Executor
import biskit.tools as t

class Ranch( Executor ):
   

    chains = {}
    pdbs_in = ()
    sequence = ''

"""
    def __init__(self, *domains, chains={}):

        # Create temporary folder for pdbs and sequence
        tempdir = tempfile.mktemp( '', self.__class__.__name__.lower() + '_', t.tempDir() )

        # Create temporary folder for models
        self.dir_models = tempfile.mktemp( '', 'models_', tempdir )

        # Temporary file for sequence
        self.f_seq = tempfile.mktemp('_sequence.seq', '', tempdir)
        
        # Generate by default 10 models, with no intensities
        args = self.f_seq + ' -q=10 -i'

        args = args + ' -w=%s' % (self.tempmodels)

        self.domains = domains
        self.sequence = ''

        Executor.__init__(self, 'ranch', tempdir=tempdir, args=args, cwd=tempdir)

    def prepare(self):
        """
        Overrides Executor method.
        """
        # Create tempdir for pdbs and seq
        Executor.prepare(self)

        # Create temporary directory for models
        if not os.path.exists(self.tempmodels):
            os.mkdir( self.tempmodels )

        ######## DIAGRAM STARTS ########

        for element in self.domains:
            if isinstance(element, str):
                self.sequence += element
            else:
                if element.symcore:   ### CREATE THIS VARIABLE... FIND A WAY TO MAKE THIS CHECK ONLY ONCE??
                    self.sequence += element
                    self.sym = 'p' + string(self.countchains())



    def cleanup(self):
        """
        Delete temporary files
        """
        if not self.debug:
            t.tryRemove(self.f_seq)
            for file in self.pdbs_in:
                t.tryRemove(file)
        Executor.cleanup(self)


    def finish(self):
        """
        Overrides Executor method.
        """
        Executor.finish(self)

        ### Retrieve models created as PDBModels
        m_paths = [os.path.join(self.tempmodels, f) for f in os.listdir(self.tempmodels)]
        self.result = [B.PDBModel(m) for m in m_paths]

      ## Fix the models
"""

#############
##  TESTING        
#############
import biskit.test as BT

class TestRanch(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]

    ## Write tests for each case in the ranch examples folder

    def prepare(self):
        # Create PDBModels from pdb files
        dom1 = B.PDBModel("/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/2z6o_mod.pdb")
        dom2 = B.PDBModel("/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/Histone_H3.pdb")
        domAB1 = B.PDBModel("/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
        domAB2 = domAB1.clone()
        
        # Call ranch for exmaple 1 in ranch_examples
        call_2z60_linker_histone = r.Ranch(dom1, 'GGGGGGGGGG', dom2)
        models1 = call_2z60_linker_histone.run()
        
        # Call ranch for example 4 in ranch_examples
        call_domAB1_linker_domAB2 = r.Ranch(domAB1, 'GGGGGGGGGG', domAB2, 
            chains = {domAB1:'A', domAB2:'B'})
        models2 = call_domAB1_linker_domAB2.run()

        # example 5 in ranch_examples
        call_domAB1_linker_domAB2 = r.Ranch(domAB1, 'GGGGGGGGGG', domAB2, 
            chains = {domAB1:'A', domAB2:'B'}, symmetry=)
        models2 = call_domAB1_linker_domAB2.run()

    def test_input_files_created(self):



    def test_pdbmodels_saved(self):
        # models lists contain 10 elements
        self.assertTrue(len(models1)==10, "models1 does not contain 10 elements")
        self.assertTrue(len(models1)==10, "models2 does not contain 10 elements")
        
        # elements of models lists are PDBModels
        self.assertTrue(isinstance(models1[0], B.PDBModel), 
            "models1 contents are not PDBModels")
        self.assertTrue(isinstance(models2[0], B.PDBModel), 
            "models2 contents are not PDBModels")

    def test_number_of_chains(self):


    def test_order_amino_acids(self):



if __name__ == '__main__':

    BT.localTest(debug=False)
