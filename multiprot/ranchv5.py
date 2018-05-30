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


    def __init__(self, **kwargs):

        # Create temporary folder for pdbs and sequence
        tempdir = tempfile.mktemp( '', self.__class__.__name__.lower() + '_', t.tempDir() )

        # Create temporary folder for models
        self.tempmodels = tempfile.mktemp( '', 'models_', tempdir )

        # Temporary file for sequence
        self.f_seq = tempfile.mktemp('_sequence.seq', '', tempdir)
        
        # Generate by default 10 models, with no intensities
        args = self.f_seq + ' -q=10 -i'

        args = args + ' -w=%s' % (self.tempmodels)

        Executor.__init__(self, 'ranch', tempdir=tempdir, args=args, cwd=tempdir)

    
    def addChain(self, id, chain):

        '''
        Add chain to the model to be generated
        @ param id: id of the chain
        @ type id: str
        @ param chain: sequence of the elements constituting the chain
        @ type chain: list or tuple (domain1, linker1, domain2...)
        @ param fixed: is each domain to remain fixed in the original
            coordinates?
        @ type fixed: list ir tuple with values ('yes', 'no') for each domain
            default: None, only first domain is fixed
        '''

        self.chains[id:(chain,fixed)]

        for elem in chain:
            if isinstance(elem, B.PDBModel):
                # Make temp pdb file
                self.pdbs_in = self.pdbs_in + (tempfile.mktemp(
                    elem.sourceFile()[-8:], '', tempdir),)
                # Add sequence to seq
                # TEST FOR NON AA CHARACTERS
                self.sequence += elem.sequence()
                
            elif isinstance(elem, str):
                self.sequence += elem
            else:
                pass

        args = args + ' -x=%s' * len(self.pdbs_in) % self.pdbs_in

        if fixed:
            args = args + ' -f=%s' * len(fixed) % fixed
        else:
            args = args + ' -f=%s' * len(self.pdbs_in) % 
            ('yes',)+('no',)*(len(self.pdbs_in)-1)  # only the first is 'yes'

    def symmetry(self, sym, only_symmetric=None):
        self.symmetry = sym
        # MODIFY TO SET SYMMETRY MORE THAN ONCE
        args = args + ' -s=%s' % self.symmetry

        if only_symmetric:
            args = args + ' -y=S'


    def multichain(self, groups):
        '''
        Define the multichain domains
        @ param groups: groups of two or more chains that will remain together
            in space
        @ type groups: tuple of tuples with two or more domain names
        '''
        self.groups = groups

    def prepare(self):
        """
        Overrides Executor method.
        """
        # Create tempdir for pdbs and seq
        Executor.prepare(self)

        # Create temporary directory for models
        if not os.path.exists(self.tempmodels):
            os.mkdir( self.tempmodels )

        with open(self.f_seq, 'w') as f:
            f.write(self.sequence)

        if self.x:
            for i in range(len(self.x)):
                self.x[i].writePdb(self.pdbs_in[i])

    def cleanup(self):
        """
        Delete temporary files
        """
        Executor.cleanup(self)
        if not self.debug:
            t.tryRemove(self.f_seq)
            for file in self.pdbs_in:
                t.tryRemove(file)

    def finish(self):
        """
        Overrides Executor method.
        """
        Executor.finish(self)

        ### Retrieve models created as PDBModels
        m_paths = [os.path.join(self.tempmodels, f) for f in os.listdir(self.tempmodels)]
        self.result = [B.PDBModel(m) for m in m_paths]

        ## Fix the models
        