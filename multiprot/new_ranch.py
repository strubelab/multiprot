"""
Author: Francisco Javier Guzman

Last edited: 20/feb/18

"""


import biskit as B
import numpy as N
import tempfile, os
import re

from biskit.exe.executor import Executor
import biskit.tools as t

class Ranch( Executor ):
   

    chains = {}
    pdbs_in = ()
    sequence = ''


    def __init__(self, *domains, chains={}, symmetry='p1', symtemplate=None, fixed=None, multich=None):
        # Raise exception if symmetry is different than p1 and symtemplate is None,
        # and if symtemplate is not in domains
        
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
        self.symmetry = symmetry
        self.symtemplate = symtemplate
        self.doms_in = []
        self.fixed = fixed
        self.multich = multich
        self.embedded = []

        if self.fixed == None:
            self.fixed = []
            for element in self.domains:
                if isinstance(element, B.PDBModel):
                    self.fixed.append('no')

        if self.multich == None:
            self.multich = []
            for element in self.domains:
                if isinstance(element, B.PDBModel):
                    if element == self.symtemplate:
                        self.multich.append('yes')
                    else:
                        self.multich.append('no')

        Executor.__init__(self, 'ranch', tempdir=tempdir, args=args, cwd=tempdir)


    def _setup(self):
        """
        - Creates the sequence from the domains and linkers
        - Creates new PDBModels for ranch input if necessary (which are later
           converted into pdb files by prepare method)
        """
        
        ######## DIAGRAM STARTS ########

        i = 0   # Counter for domain position
        for element in self.domains:
            if isinstance(element, str):    # 1
                # If is sequence, add to sequence and continue
                self.sequence += element
                continue

            elif isinstance(element, B.PDBModel):
                # if is PDBModel

                if self.multich[i]=='yes' or element.lenChains()==1:    # 2
                    # If is symmetry core or single domain
                    # Find a way to make the symmetry test only once?
                    # Action: Add to sequence and to domains list
                    self.sequence += element.sequence()
                    self.doms_in.append(element)
                    i += 1

                else:  
                    # Not symmetric, and part of a multichain complex
                    ## Look only at domain that corresponds to this chain

                    # Get chain index
                    chain_id = chains[element]
                    mask_chain = element.maskFrom('chain_id', chain_id)
                    i_mask_chain = N.nonzero(mask_chain)[0]
                    chain_ind = element.atom2chainIndices(i_mask_chain)[0]

                    if self.fixed[i] == 'yes':  # 3
                        # Already modeled, paired with another chain
                        # THE HIGHER LEVEL PROGRAM MUST GIVE THE ENTIRE MODEL
                        # AS ONE OF THE FIXED DOMAINS, TO AVOID STEPS 4 AND 5
                        # IN DIAGRAM
                        # CONTINUE IN NEW_RANCH_2 WITHOUT STEP 5

                        if element in self.embedded:    # 5
                            # If the domain is included in the embedded domains,
                            # meaning that the paired chain is already embedded
                            # in a previous domain ... is this the best way?

                            # Action: Take fixed domain only as a separate
                            # PDBModel, add to sequence and domains list
                            m = element.takeChains([chain_ind])
                            self.sequence += m.sequence()
                            self.doms_in.append(m)

                        else:
                            # The paired chain is not embedded in a previous domain
                            m = element.takeChains([chain_ind])
                            embedded = self.embed(m,element)    

                            # AT THE END SAVE HOW THE SEQUENCE SHOULD LOOK, TO HELP FIND SYMMETRIC UNIT

            else:
                raise TypeError(
                    'The *domains arguments must be either strings or PDBModels')

    @staticmethod
    def embed(self, dom, full):
        '''
        Embeds one model (full - possibly with multiple chains) into another 
        (dom) to trick ranch to treat them as a simple single-chain domain

        :param dom: model of a single chain domain that will contain the
                       other model
        :type dom: PDBModel
        :param full: model of a single or multiple chain domain, that generally
                     contains also the 'dom' domain, and will be embedded
                     into 'dom'
        :type full: PDBModel
        '''

        # Take dom from full
        start = dom.sequence()[0:10]
        if re.search(start, full.sequence()):
            # If the dom sequence is inside full sequence
            matches = re.finditer(start, full.sequence())
            first_res_dom = dom.res2atomindices([0])
            lowdom = first_res_dom[0]
            highdom = first_res_dom[-1]
            
            for match in matches:
                index = match.start()
                first_res_full = full.res2atomindices([index])

                lowfull = first_res_full[0]
                highfull = first_res_full[-1]

                if N.all(dom.xyz[lowdom:highdom+1] == full.xyz[lowfull:highfull+1]):
                    # If the atoms for the first residue are in the same positions
                    # Take chain
                    chain_ind = full.atom2chainIndices(first_res_full)
                    chains_to_take = list(range(full.lenChains()))
                    chains_to_take = chains_to_take.remove(chain_ind[0])

                    full = full.takeChains(chains_to_take)


    def prepare(self):
        """
        Overrides Executor method.
        """
        # Create tempdir for pdbs and seq
        Executor.prepare(self)

        # Create temporary directory for models
        if not os.path.exists(self.tempmodels):
            os.mkdir( self.tempmodels )



    def cleanup(self):
        """
        Delete temporary files
        """
        if not self.debug:
            t.tryRemove(self.f_seq)
            for file in self.pdbs_in:
                t.tryRemove(file)

        super().cleanup(self)

    def finish(self):
        """
        Overrides Executor method.
        Write more....
        """
        # CHECK IF IS CALLED ONLY IF SUCCESSFUL
        Executor.finish(self)

        ### Retrieve models created as PDBModels
        m_paths = [os.path.join(self.tempmodels, f) for f in os.listdir(self.tempmodels)]
        self.result = [B.PDBModel(m) for m in m_paths]


      ## Fix the models


#############
##  TESTING        
#############
import biskit.test as BT

class TestRanch(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE, BT.LONG ]

    ## Write tests for each case in the ranch examples folder

    def prepare(self):
        pass

    def test_input_files_created(self):
        pass


    def test_pdbmodels_saved(self):

        # Create PDBModels from pdb files
        dom1 = B.PDBModel(
            "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/2z6o_mod.pdb")
        
        dom2 = B.PDBModel(
            "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/Histone_H3.pdb")
        
        domAB1 = B.PDBModel(
            "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
        
        domAB2 = domAB1.clone()
        
        # Call ranch for exmaple 1 in ranch_examples
        call_example1 = Ranch(dom1, 'GGGGGGGGGG', dom2)
        models1 = call_example1.run()
        
        # Call ranch for example 4 in ranch_examples
        call_example4 = Ranch(domAB1, 'GGGGGGGGGG', domAB2, 
            chains = {domAB1:'A', domAB2:'B'})
        models2 = call_example4.run()

        # example 5 in ranch_examples
        call_example5 = Ranch(domAB1, 'GGGGGGGGGG', domAB2, 
            chains = {domAB1:'A', domAB2:'B'}, symmetry='p2', symtemplate=domAB1)
        models3 = call_example5.run()

        # models lists contain 10 elements
        self.assertTrue(len(models1)==10, "models1 does not contain 10 elements")
        self.assertTrue(len(models2)==10, "models2 does not contain 10 elements")
        self.assertTrue(len(models3)==10, "models3 does not contain 10 elements")
        
        # elements of models lists are PDBModels
        self.assertTrue(isinstance(models1[0], B.PDBModel), 
            "models1 contents are not PDBModels")
        self.assertTrue(isinstance(models2[0], B.PDBModel), 
            "models2 contents are not PDBModels")
        self.assertTrue(isinstance(models3[0], B.PDBModel), 
            "models1 contents are not PDBModels")



class TestCleaning(BT.BiskitTest):
    """ Test class for the cleaning methods post-run """

    def test_number_of_chains(self):
        pass


    def test_order_amino_acids(self):
        pass



if __name__ == '__main__':

    BT.localTest(debug=False)
