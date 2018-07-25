#RG: please convert to 4-space indentation (e.g. WingIDE can do it automatically)
#RG: remove non-ascii characters otherwise we need to specify encoding
"""
Author: Francisco Javier Guzman-Vega

Wrapper for RANCH (RANdom CHain) - tool for the generation of a pool of random
models based upon user supplied sequence and structural information

https://www.embl-hamburg.de/biosaxs/eom.html
https://www.embl-hamburg.de/biosaxs/manuals/eom.html

Tria, G., Mertens, H. D. T., Kachala, M. & Svergun, D. I. (2015) Advanced 
    ensemble modelling of flexible macromolecules using X-ray solution scattering. 
    IUCrJ 2, 207-217

"""

## STILL NEED TO WRITE CODE TO TEST AND HANDLE INCORRECT INPUT, AND WHAT TO DO
## IF RANCH FAILS

## DELETE SELF. VARIABLES? (E.G. SELF.DOMS_IN)

import biskit as B
import numpy as N
import tempfile, os
import re
from operator import itemgetter
from errors import *
from biskit.exe.executor import Executor
import biskit.tools as T #RG: I prefer uppercase (T) to distinguish variables and modules (SOLVED)

#### Helper tool misc functions ####

def embed(dom, to_embed):
    """
    Embeds one model (int_dom - possibly with multiple chains) into another 
    (dom) to trick ranch to treat them as a simple single-chain domain

    :param dom: model of a single chain domain that will contain the
                    other model
    :type dom: PDBModel
    :param int_dom: model of a single or multiple chain domain that will be 
                     embedded into 'dom'
    :type int_dom: PDBModel
    """

    first = dom.take(dom.res2atomIndices([0,1]))
    last = dom.take(dom.res2atomIndices(list(range(2,dom.lenResidues()))))

    return first.concat(to_embed,last)


def extract_fixed(dom, full):
    """
    Extracts one model from another
    Finds the position of 'dom' inside 'full' comparing the sequence and atom
    coordinates for each chain in dom, gets the chain index and takes all the
    chains but the ones selected.
    
    :param dom: model of a single or multiple chain domain
    :type dom: PDBModel
    :param full: model of a multiple chain domain that contains 'dom'
    :type full: PDBModel

    :return: model 'full' without dom
    :type return: PDBModel
    """

    chains_to_take = list(range(full.lenChains()))

    # Make a list with one PDBModel for each chain in dom
    # This is to find one chain from dom at a time, in case they are not
    # together in 'full' ... is this even necessary?
    doms = [dom.takeChains([i]) for i in range(dom.lenChains())]

    for m in doms:

        start = m.sequence()[:10]  # could use the entire sequence instead
        
        if re.search(start, full.sequence()):
            # If the m sequence is inside full sequence
            # Action: look for the position of m inside full, and extract

            matches = re.finditer(start, full.sequence())
            first_res_m = m.res2atomIndices([0])
            lowm = first_res_m[0]
            highm = first_res_m[-1]
            
            for match in matches:
                index = match.start()
                first_res_full = full.res2atomIndices([index])

                lowfull = first_res_full[0]
                highfull = first_res_full[-1]

                if N.all(m.xyz[lowm:highm+1] == full.xyz[lowfull:highfull+1]):
                    # If the atoms for the first residue are in the same positions
                    # Action: remove chain index from chains_to_take
                    chain_ind = full.atom2chainIndices(first_res_full)
                    chains_to_take.remove(chain_ind[0])
                    break

    full = full.takeChains(chains_to_take)

    return full


def extract_embedded(full, embedded):
    """
    Extracts one  or more PDBModels from another
    Finds the sequence and location of each domain in embedded dictionary.
    Extracts the atoms and concatenates at the end of 'self'. 
    Renumbers amino acids, id number and renames chains in the process.
    
    :param embedded: dictionary with embedded domains and its position (index)
                            in the full sequence
    :type embedded: dictionary

    :return: 'full' with embedded domains concatenated at the end as 
                independent chains
    :type return: PDBModel
    """

    ## For the higher level program, add an argument to provide the dictionary

    chains_to_take = list(range(full.lenChains()))

    r = B.PDBModel()
    emb_ind = []   # List for start and end indexes for each embedded domain

    for dom, i_start in embedded.items():
        
        i_end = i_start + len(dom.sequence())

        if full.sequence()[i_start:i_end] == dom.sequence():
            r = r.concat(full.takeResidues(list(range(i_start, i_end))))
            emb_ind.append((i_start, i_end))
        else:
            raise MatchError('The sequence from the domain to exctract does not \
                match the sequence in the full domain with the specified indices')

    # Sort the list to remove the atoms from highest to lowest index,
    # so the indexes won't be affected
    emb_ind = sorted(emb_ind, key=itemgetter(0), reverse=True)

    for i_start, i_end in emb_ind:
        atomi_start = full.resIndex()[i_start]
        atomi_end = full.resIndex()[i_end]
        full.remove(list(range(atomi_start, atomi_end)))

    # Concat the original chain that previously contained the embedded domains
    for i in range(full.lenChains()-1):
        full.mergeChains(0)

    full.renumberResidues()     # Renumber amino acids
    full = full.concat(r)       # Combine full and r
    full.addChainId()           # Add chain IDs with consecutive letters
                                # NOTE: add feature for personalized chain names

    # Renumber atoms
    full['serial_number'] = N.arange(1,len(full)+1)

    return full


def extract_symmetric(full, symseq, embedded):
    """
    Extracts one or more embedded chains from a PDBModel with a symmetric
    structure
    
    :param full:   PDBModel with symmetric structure, that contains embedded 
                        chains
    :type full: PDBModel

    :param symseq: sequence of the symmetric unit, i.e. the sequence that is
                        multiplied in the symmetric structure
    :type symseq: string
    :param embedded: dictionary with embedded domains and its position (index)
                            in the sequence of 'full'
    :type embedded: dictionary
    :return: 'full' with embedded domains concatenated at the end for each
                symmetric unit
    :type full: PDBModel
    """
    
    symunits = []

    if re.search(symseq, full.sequence()):

        matches = re.finditer(symseq, full.sequence())

        for match in matches:
            istart, iend = match.span()
            symunit = full.takeResidues(list(range(istart, iend)))
            # Extract embedded domains one symunit at a time
            symunits.append(extract_embedded(symunit, embedded))

        r = symunits[0]

        for i in range(1,len(symunits)):
            r = r.concat(symunits[i])

        r.addChainId()
        r['serial_number'] = N.arange(1,len(r)+1)

    else:
        raise MatchError("Symseq could not be found inside the full domain")

    return r


class Ranch(Executor):

    """
    A Ranch wrapper to generate 10 independent models based on sequence and
    protein structures in the form of PDBModels.

    Usage
    =====

    >>> call = Ranch(dom1, linker1, dom2, linker2, ..., chains = {dom1:'chainA',
                        dom2:'chainC'}, symmetry='p2', symtemplate=dom2, 
                        pool_sym='symmetry')
    >>> models = call.run()

    The wrapper takes the sequence of domain and linker elements that will make
    up the models and returns a list of 10 models as PDBModel objects.

    """


    def __init__(self, *domains, chains={}, symmetry='p1', symtemplate=None, 
        symunit=None, pool_sym='m', fixed=[], n=10, **kw):
        
        """
        Creates the variables that Ranch needs to run
        
        ## There is no need to provide chains dict entry for symtemplate

        :param *domains:  Domains and linkers that will make up the chain to be 
                                modeled, in succession
        :type *domains:   PDBModels for structured domains and strings for the
                                linkers, separated by commas
        :param chains:  Specifies which chain will be taken from each domain if
                        they have multiple chains. If this parameter is ommited
                        the program will take the first chain by default
        :type chains:   Dictionary with entries in the form of {domain:chain_id}
                        where 'domain' is a PDBModel and 'chain_id' is a string
        :param symmetry: Specifies the type of symmetry that the model will have.
                         E.g. 'p1' is for no symmetry, 'p2' for 'duplicate' 
                         symemtry (dimers), 'p3' for 'triplicate' symmetry 
                         (tetramers), etc.
        :type symmetry:   string
        :param symtemplate:  Specifies the domain that will be taken as symmetry
                             core, in case of a symmetric model being created.
                             The symtemplate has to be a symmetric molecule 
                             itself. The number of chains in the symtemplate has 
                             to agree with the type of symmetry chosen. E.g., if 
                             symmetry is 'p2', symtemplate has to be a dimer. If 
                             symmetry is 'p3', symtemplate should be a trimer, 
                             etc.
        :type symtemplate:   PDBModel
        :param symunit: Sequence of the symmetric unit (sequence that is 
                        repeated) in the symtemplate, if it consists of multiple
                        chains. This parameter will probably be used only by the
                        higher level implementation of the wrapper.
        :type symunit:  string
        :param pool_sym:  Specifies the overall symmetry of the models 
                          generated. 'symmetry' will generate symmetric 
                          multichain structures only. 'asymmetry' will generate
                          models with a symmetric core structure but leaves the
                          remaining structure asymmetric, and 'mix' will create
                          both symmetric and asymmetric models. This argument 
                          is not used when the symmetry is 'p1' (no symmetry). 
                          If symmetry is different to 'p1' and pool_sym is 
                          not specified it defaults to 'mix'.
        :type pool_sym:   string
        :param fixed: Specifies if the original coordinates are to be maintained
                      for one or more domains (PDBModels) in the input. This 
                      argument will only be used by the higher level 
                      implementation of the wrapper, except in special cases when
                      the user wishes to maintain the coordinates of a specific 
                      domain. If more than one domain is specified as fixed, it 
                      may cause an error if the linker length is not appropriate 
                      to bind these domains. If no domain is specified as fixed, 
                      ranch automatically fixes the symmetry core if present; if 
                      not, it fixes the first domain only.
        :type fixed: List with the domains (PDBModels) to be maintained as fixed
        :param n: Number of models to be generated (let's say from 10 to 15,000,
                  though I'm not sure about the actual maximum for ranch)
        :type n:    Integer
        :param kw:  additional key=value parameters are passed on to
                    'Executor.__init__'. For example:
                    ::
                        debug    -  0|1, keep all temporary files (default: 0)
                        verbose  -  0|1, print progress messages to log
                                        (log != STDOUT)
                        nice     -  int, nice level (default: 0)
                        log      -  biskit.LogFile, program log (None->STOUT)
                                        (default:None)

        """

        # TODO: Raise exception if symmetry is not p1 and symtemplate is 
        # None, and if symtemplate is not in domains

        #TODO: Add possibility to input more options for ranch
        
        # Create temporary folder for pdbs and sequence
        #RG: Executor can do that for you if you set `tempdir` parameter to True or to a custom name
        #JG:  I needed to create this before calling Executor.__init__(...) to make 
        #     another folder inside (line 318)
        tempdir = tempfile.mkdtemp( '', self.__class__.__name__.lower() + '_', 
            T.tempDir() )

        # Create temporary folder for models
        self.dir_models = tempfile.mkdtemp( '', 'models_', tempdir )

        #RG: preferred (cross-platform): os.path.join(tempdir, 'sequence.seq') (SOLVED)
        self.f_seq = os.path.join(tempdir, 'sequence.seq')

        self.n = n
        self.domains = domains
        self.chains = chains
        self.sequence = ''
        self.symmetry = symmetry
        self.symtemplate = symtemplate
        self.doms_in = []    # list of domains as PDBModels
        self.pdbs_in = []    # list of pdb file paths
        self.embedded = {}   # dictionary with domain : residue number to
                                    # identify and locate embedded domains

        self.pool_sym = pool_sym

        #RG: list comprehension should work instead of the complex loop (not tested)
        #JG:  good idea, modified it a bit because the list still needs to contain
        #     'yes' or 'no' values to be pasted into the args string for ranch
        
        # self.fixed = [ element in fixed for element in self.domains if \
        #               isinstance(element, B.PDBModel) ]

        self.fixed = ['yes' if element in fixed else 'no' for element in \
                            self.domains if isinstance(element, B.PDBModel)]

        self.multich = ['yes' if element is self.symtemplate else 'no' for \
                                element in self.domains if isinstance(element, 
                                    B.PDBModel)]

        if symtemplate:
            if symunit:
                self.symunit = symunit
            else:
                # If there is no symunit provided, it is a single chain symunit
                # Action: take symunit from symtemplate
                self.symunit = symtemplate.takeChains([0])

        # Path for config file
        self.configpath = [os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'exeConfig/')]

        super().__init__('ranch', tempdir=tempdir, cwd=tempdir,
            configpath=self.configpath, **kw)

    def _setup(self):
        """
        - Creates the sequence from the domains and linkers
        - Creates new PDBModels for ranch input if necessary (which are later
            converted into pdb files by prepare method)
        """
        
        ######## DIAGRAM STARTS ########

        ## NOTE: The numbers are references to the steps in DIAGRAM.png

        i = 0   # Counter for domain position. Counts only PDBModels
        for element in self.domains:
            
            if isinstance(element, str):    # 1
                # If is sequence, add to sequence and continue
                self.sequence += element
                continue      #RG: continue and break are, if at all possible, not supposed to be used in good code :)

            elif isinstance(element, B.PDBModel):
                # if is PDBModel

                if self.multich[i]=='yes':
                    # If is symtemplate, there is symmetry ... 2

                    # For multiple chains with symmetry, the symetric unit should 
                    # already be embedded in the (modified) symtemplate, and the
                    # symmetric unit by itself provided as the argument symunit

                    # Action: Add sequence of symunit to sequence and the element
                    # to domains list

                    # Find a way to make the symmetry test only once?

                    self.sequence += self.symunit.sequence()
                    self.doms_in.append(element)

                elif element.lenChains()==1 or self.fixed[i] == 'yes':
                    # If is single chain-domain ... 2.1
                    # Or is a domain already modeled/fixed ... 3

                    # THE HIGHER LEVEL PROGRAM MUST GIVE THE ENTIRE MODEL WITH
                    # MULTIPLE CHAINS EMBEDDED IN ONE OF THE FIXED DOMAINS, AND
                    # IT WILL BE TREATED AS A SINGLE CHAIN-DOMAIN, TO AVOID STEPS 
                    # 4 AND 5 IN DIAGRAM
                    # IT ALSO HAS TO REMOVE THE BOUND DOMAINS THAT WILL BE
                    # MODELED IN OTHER CHAINS, TO AVOID STEP 6
                    
                    # Action: Add to sequence and to domains list
                    # Conserve fixing
                    
                    self.sequence += element.sequence()
                    self.doms_in.append(element)

                else:
                    # Not modeled, but paired with another domain

                    if element in self.chains and \
                        isinstance(self.chains[element], (list, tuple)):   # 7
                        # if the other part of the domain is bound to the
                        # same chain, the chains entry of the domain will
                        # have a list or tuple with the chains to be taken, 
                        # in order

                        # Action: Embed sequence of paired domain into
                        # element. Take the sequence and domains up to this
                        # point and model. Change coordinates, fix both
                        # domains and separate as individual fixed domains. 
                        
                        # Go back before adding the first dom (or start again),
                        # proceed as normal
                        # ALTERNATIVELY, take modeled region up to first dom as 
                        # fixed?

                        # Be careful to put the separate domains in the correct 
                        # order, consulting self.chains[element]

                        # Can only work in the higher level implementation

                        # 8
                        # Get chain index
                        chain_id = self.chains[element][0]
                        # Gets True values only for the specified chain index
                        mask_chain = element.maskFrom('chain_id', chain_id)
                        # Convert True values to indices
                        i_mask_chain = N.nonzero(mask_chain)[0]
                        chain_ind = element.atom2chainIndices(i_mask_chain)[0]

                        m = element.takeChains([chain_ind])
                        to_embed = extract_fixed(m, element)
                        m_emb = embed(m, to_embed)

                        self.embedded[to_embed] = len(self.sequence) + 2

                        self.sequence += m_emb.sequence()
                        self.doms_in.append(m_emb)

                        break

                    else: #RG: mhm... what is this embedding thing about? ... MAGIC
                        # Only one domain from element is part of the chain
                        # Action: Embed the paired domains into the selected chain

                        # 9
                        # Get chain index
                        if element in self.chains:
                            chain_id = self.chains[element]
                            mask_chain = element.maskFrom('chain_id', chain_id)
                            i_mask_chain = N.nonzero(mask_chain)[0]
                            chain_ind = element.atom2chainIndices(i_mask_chain)[0]
                        else:
                            chain_ind = 0

                        m = element.takeChains([chain_ind])
                        to_embed = extract_fixed(m, element)
                        m_emb = embed(m, to_embed)
                        
                        self.embedded[to_embed] = len(self.sequence) + 2

                        self.sequence += m_emb.sequence()
                        self.doms_in.append(m_emb)

        ####### DIAGRAM FINISHES... REACHED STEP 9 #########

            else:  
                #RG: I think it's considered better to raise your own custom errors rather than python built-in (SOLVED)
                raise InputError(
                    'The *domains arguments must be either strings or PDBModels.')

            i += 1

        # Symseq is the sequence that will be multiplied in the symmetric 
        # structure
        if self.symtemplate:
            self.symseq = self.sequence

        return None


    def prepare(self):
        """
        Overrides Executor method.
        """

        try:
            self._setup()
        except:
            self.cleanup()
            raise

        # Write pdb files
        for i in range(len(self.doms_in)):
            pdb_name = os.path.join(self.tempdir, str(i)+'_')
            if self.doms_in[i].validSource() is None:     #RG: slightly faster: ...validSource() is None (SOLVED)
                # If it was a pdb created 'de novo'
                pdb_name += '.pdb'
            else:
                # If it comes directly from a file
                pdb_name += self.doms_in[i].sourceFile()[-8:]

            self.pdbs_in.append(pdb_name)
            self.doms_in[i].writePdb(pdb_name)

        # Write sequence file
        with open(self.f_seq, 'w') as f:
            f.write(self.sequence)

        # Generate n models with no intensities
        self.args = self.f_seq + ' -q=%s -i' % self.n   #RG: make this another __init__ parameter! (SOLVED)

        if self.symtemplate:
            self.args = self.args + ' -s=%s -y=%s' % (self.symmetry, 
                self.pool_sym)

        self.args = self.args + ' -x=%s' * len(self.pdbs_in) % tuple(self.pdbs_in)

        self.args = self.args + ' -f=%s' * len(self.fixed) % tuple(self.fixed)

        self.args = self.args + ' -o=%s' * len(self.multich) % tuple(self.multich)

        self.args = self.args + ' -w=%s' % self.dir_models

        
    def isFailed(self):
        """
        Overrides executor method
        """
        ## CHECK HOW TO SEE OUTPUT

        return self.error or 'Problems' in str(self.output)

    def fail(self):
        """
        Overrides Executor method. Called when execution fails.
        """

        ## CHECK IF THE MESSAGE POPS UP

        print('Ranch call failed.')   ## temporary

        ## PRINT ERROR MESSAGE FROM RANCH


    def finish(self):
        """
        Overrides Executor method.
        Write more....
        """
        
        ### Retrieve models created as PDBModels
        m_paths = [os.path.join(self.dir_models, f) for f in os.listdir(
            self.dir_models)]

        if self.symtemplate:
            self.result = [extract_symmetric(B.PDBModel(m), self.symseq,
                self.embedded) for m in m_paths]
        else:
            self.result = [extract_embedded(B.PDBModel(m), self.embedded
                ) for m in m_paths]


    def cleanup(self):
        """
        Delete temporary files
        """
        if not self.debug:
            T.tryRemove(self.tempdir, tree=True)




#############
##  TESTING        
#############
import testing

class TestRanch(testing.AutoTest):
    """
    Test class
    Test examples 1, 4, 5, 7, and 10 from ranch_examples/
    
    Runs each case separately and tests for the creation of the models list
    with 10 PDBModels, and then for one of the models tests for number of chains,
    length of each chain, residue numbering and serial numbering

    These tests take around 35 seconds to complete, because ranch is called for
    each example
    """

    TAGS = [ testing.EXE, testing.LONG ]

    #RG: I know this is convenient but very bad idea to execute any code in the class definition body (SOLVED)
    
    #RG: problem 1: this will need to be executed whenever the module is loaded (not just for testing)
    #JG:  - why run it every time the module is loaded?
    #     - these tests take like 37 seconds
    
    #RG: problem 2: this will break on any other computer except your own (paths :) ) (SOLVED(?))

    #RG: one pattern that should work instead:
    dom1 = None ## define empty class variable
    dom2 = None 
    domAB1 = None
    domAB2 = None
    testpath = None

    def setUp(self):
        # self.DOM1 = DOM1 or B.PDBModel( T.testRoot('ranch/1/2z6o_mod.pdb') )
        # self.DOM2 = DOM2 or B.PDBModel( T.testRoot('ranch/1/Histone_H3.pdb') )
        ## this will load the PDBs only once even though setup is run for every test
        ## doesn't it need the 'self.'DOM1 ?
        
        ## Is this enough so it doesn't break in other computers?
        self.testpath = self.testpath or \
            os.path.join(os.path.abspath(os.path.dirname(__file__)), 'testdata')
        
        self.dom1 = self.dom1 or B.PDBModel(os.path.join(self.testpath, 
            '2z6o_mod.pdb'))
        self.dom2 = self.dom2 or B.PDBModel(os.path.join(self.testpath, 
            'Histone_H3.pdb'))
        self.domAB1 = self.domAB1 or B.PDBModel(os.path.join(self.testpath, 
            'dom1_AB.pdb'))
        self.domAB2 = self.domAB2 or self.domAB1.clone()

    def test_example1(self):
        call = Ranch(self.dom1,'GGGGGGGGGG',self.dom2)
        models = call.run()

        self.assertTrue(len(models)==10, "models does not contain 10 elements")
        self.assertTrue(isinstance(models[0], B.PDBModel), 
            "models contents are not PDBModels")

        model = models[0]
        self.assertTrue(model.lenChains()==1, 'Incorrect number of chains')
        self.assertTrue(len(model.sequence())==274, 'Incorrect chain length')
        self.assertTrue(model.atoms['residue_number'][-1]==274, 
            'Incorrect residue numbering')
        self.assertTrue(model.atoms['serial_number'][-1]==2181, 
            'Incorrect serial numbering')

    def test_example4(self):
        call = Ranch(self.domAB1, 'GGGGGGGGGGGGGGGGGGGG', self.domAB2, 
            chains = {self.domAB1:'A', self.domAB2: 'B'})
        models = call.run()

        self.assertTrue(len(models)==10, "models does not contain 10 elements")
        self.assertTrue(isinstance(models[0], B.PDBModel), 
            "models contents are not PDBModels")

        model = models[0]
        self.assertTrue(model.lenChains()==3, 'Incorrect number of chains')
        self.assertTrue(len(model.takeChains([0]).sequence())==456 and \
            len(model.takeChains([1]).sequence())==218 and \
            len(model.takeChains([2]).sequence())==218, 'Incorrect chain length')
        self.assertTrue(model.takeChains([0]).atoms['residue_number'][-1]==456 \
            and model.takeChains([1]).atoms['residue_number'][-1]==218 \
            and model.takeChains([2]).atoms['residue_number'][-1]==218, 
            'Incorrect residue numbering')
        self.assertTrue(model.atoms['serial_number'][-1]==7163, 
            'Incorrect serial numbering')

    def test_example5(self):
        call = Ranch(self.domAB1, 'GGGGGGGGGGGGGGGGGGGG', self.domAB2, 
            chains = {self.domAB2: 'A'}, symmetry='p2', symtemplate=self.domAB1, 
            pool_sym='s')
        models = call.run()

        self.assertTrue(len(models)==10, "models does not contain 10 elements")
        self.assertTrue(isinstance(models[0], B.PDBModel), 
            "models contents are not PDBModels")

        model = models[0]
        self.assertTrue(model.lenChains()==4, 'Incorrect number of chains')
        self.assertTrue(len(model.takeChains([0]).sequence())==456 and \
            len(model.takeChains([1]).sequence())==218 and \
            len(model.takeChains([2]).sequence())==456 and \
            len(model.takeChains([3]).sequence())==218, 'Incorrect chain length')
        self.assertTrue(model.takeChains([0]).atoms['residue_number'][-1]==456 \
            and model.takeChains([1]).atoms['residue_number'][-1]==218 \
            and model.takeChains([2]).atoms['residue_number'][-1]==456 \
            and model.takeChains([3]).atoms['residue_number'][-1]==218, 
            'Incorrect residue numbering')
        self.assertTrue(model.atoms['serial_number'][-1]==10754, 
            'Incorrect serial numbering')

    def test_example7(self):
        linker = 'GGGGGGGGGGGGGGGGGGGG'
        call = Ranch(self.dom2, linker, self.domAB1, linker, self.dom2, 
            symmetry='p2', symtemplate=self.domAB1, pool_sym='mix')
        models = call.run()

        model = models[0]
        self.assertTrue(model.lenChains()==2, 'Incorrect number of chains')
        self.assertTrue(len(model.takeChains([0]).sequence())==454 and \
            len(model.takeChains([0]).sequence())==454, 'Incorrect chain length')
        self.assertTrue(model.takeChains([0]).atoms['residue_number'][-1]==454 \
            and model.takeChains([1]).atoms['residue_number'][-1]==454, 
            'Incorrect residue numbering')
        self.assertTrue(model.atoms['serial_number'][-1]==6876, 
            'Incorrect serial numbering')

    def test_example10(self):
        call = Ranch(self.domAB1, 'GGGGGGGGGGGGGGGGGGGG', self.domAB2, 
            'GGGGGGGGGGGGGGGGGGGG', self.domAB2, chains = {self.domAB2:'B'})
        models = call.run()

        model = models[0]
        self.assertTrue(model.lenChains()==4, 'Incorrect number of chains')
        self.assertTrue(len(model.takeChains([0]).sequence())==694 and \
            len(model.takeChains([1]).sequence())==218 and \
            len(model.takeChains([2]).sequence())==218 and \
            len(model.takeChains([3]).sequence())==218, 'Incorrect chain length')
        self.assertTrue(model.takeChains([0]).atoms['residue_number'][-1]==694 \
            and model.takeChains([1]).atoms['residue_number'][-1]==218 \
            and model.takeChains([2]).atoms['residue_number'][-1]==218 \
            and model.takeChains([3]).atoms['residue_number'][-1]==218, 
            'Incorrect residue numbering')
        self.assertTrue(model.atoms['serial_number'][-1]==10754, 
            'Incorrect serial numbering')


if __name__ == '__main__':

    testing.localTest(debug=False)
