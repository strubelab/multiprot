"""
Script that contains the higher level implementation of the ranch wrapper,
to handle single and multiple chain scenarios

"""


#### NOTE: REMOVE ALL full_chain REFERENCES

##### FOR THE EARLY IMPLEMENTATION OF MULTIPROT, WE WILL ASSUME ONLY ONE CHAIN
##### IN THE INPUT

import tempfile
import os, re, operator
import biskit as B
import biskit.tools as T
import numpy as N
from operator import itemgetter
import multiprot.ranch as R
import multiprot.pulchra as P
from multiprot.errors import *


class Builder:
    """
    Class that will contain all the necessary methods to call the Ranch and
    pulchra wrappers multiple times and build single and multiple chain 
    structures. This class will be used by the multiprot script.
    
    """
    full_chains = []    # Will contain the symmetric units after each modeling step
                        # (only one symmetric unit if there is no symmetry)

    def __init__(self, chains, debug, number, dest):
        """
        :param args: Object that contains the arguments parsed from the command line
        :type args: argparse.Namespace object created by calling parser.parse_args()
        """
        self.CHAINS = chains    # Original chains and PDBModels from input

        self.debug = debug
        self.num = number
        self.dest = dest

    def find_paired(self, i):
        """
        Finds to which chain(s) is chain i paired to
        
        :param i: Index of the chain to be compared against
        :type i: int

        :return paired_to:  dictionary with indices of the chains bound to chain i
                            and the names and chain_ids of the bound pdbs:

        paired_to_i = {j:[pair_ij1,pair_ij2],k:[pair_ik1],...}

        Where 
        - i is the index of the chain whose bound chains will be found
        - j and k are indexes of bound chains
        - pair_ij1 and pair_ij2 are the names of the domains that bind chains i 
          and j, in the form:
            pair_ijx = [(pdb_namei, chain_idi),(pdb_namej, chain_idj)]
        - pair_ik1 has the names of the domain that binds chains i and k, in the
          form:
            pair_ikx = [(pdb_namei, chain_idi),(pdb_namek, chain_idk)]
        """

        paired_to = {}

        chain_idxs = list(range(len(self.CHAINS)))
        chain_idxs.remove(i)

        for j in chain_idxs:
            for pdbj, chainj in self.CHAINS[j].chains_names.items():
                for pdbi, chaini in self.CHAINS[i].chains_names.items():
                    if pdbj==pdbi and chainj!=chaini:
                        if j in paired_to:
                            paired_to[j].append([(pdbi,chaini),(pdbj,chainj)])
                        else:
                            paired_to[j] = [[(pdbi,chaini),(pdbj,chainj)]]

                    # If same pdb and chain are referenced in different chains
                    # in args['chains'] parameter, they are not bound

        return paired_to


    def extract_fixed(self, dom, full):
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

        first_res_dom = dom.res2atomIndices([0])
        lowdom = first_res_dom[0]
        highdom = first_res_dom[-1]
        i=0
        for match in re.finditer(dom.sequence(), full.sequence()):
            start,end = match.span()
            first_res_full = full.res2atomIndices([start])
            lowfull = first_res_full[0]
            highfull = first_res_full[-1]
            i += 1
            # Compare only the atoms of the first residue
            if N.all(dom.xyz[lowdom:highdom+1] == full.xyz[lowfull:highfull+1]):
                atom_start = full.resIndex()[start]
                atom_end = full.res2atomIndices([end-1])[-1] + 1
                full.remove(list(range(atom_start,atom_end)))
                break

        return full

    def chainIndex(self,model,chain_id):
        """
        Take as input the chain id of model and return its index
        """
        mask_chain = model.maskFrom('chain_id', chain_id)
        i_mask_chain = N.nonzero(mask_chain)[0]
        return model.atom2chainIndices(i_mask_chain)[0]

    def replace_modeled(self, chaini, bound_indexes, s):
        """
        MODIFY
        Replaces the domains (coordinates) already modeled in a previous
        chain in the 'domains' attribute of chain j
        
        Replaces only the domains that are bound to other chains.

        In chain j it leaves only the part of the domain belonging to that chain 
        (j_dom), and removes j_dom from the (already modeled) chain i

        It adds each j_dom to the chainj.jdomains attribute, so the symmetric
        j_doms are not lost and can be used later to embed the symmetric units

        :param chaini:  chain to be processed
        :type chaini:   multiprot.Chain object
        :param j:   index of the chain whose domains will be replaced
        :type j:    int
        :param s:   number of the symmetric units to process
        :type s:    int
        
        :return None:

        Other output/modification to external variables:
        The method will replace the corresponding domains in chainj from self.CHAINS
        """
        
        modeled_domains = chaini.modeled_domains[s]
        for j in bound_indexes:
            chainj = self.CHAINS[j]
            if not chainj.modeled:
                # For every pair of coupled domains [(namei),(namej)]
                for pair in chaini.paired_to[j]:
                    # Find the index of the domain in chains i and j
                    i_ind = chaini.names.index(pair[0])
                    j_ind = chainj.names.index(pair[1])

                    if not chainj.new_domains[j_ind] or s>0:
                        # Add domain with new coordinates to chainj.new_domains, and fix
                        # Replace with only the part belonging to chain j
                        j_id = pair[1][1]
                        modeled_dom = modeled_domains[i_ind]

                        j_chainind = self.chainIndex(modeled_dom,j_id)
                        j_dom = modeled_dom.takeChains([j_chainind])

                        if s==0:
                            chainj.new_domains[j_ind] = j_dom
                            chainj.args["fixed"].append(j_dom)

                        if j_ind in chainj.jdomains:
                            chainj.jdomains[j_ind].append(j_dom)
                        else:
                            chainj.jdomains[j_ind] = [j_dom]

        return None

    def call_ranch(self, chaini):
        """
        Builds models with ranch
        """
        # Model with ranch
        call = R.Ranch(*chaini.domains, **chaini.args, debug=self.debug)
        models = call.run()
        if models is None:
            raise RanchError('Models not produced.')
        return models


    ## FOR ONLY ONE SYMMETRIC SEQUENCE
    def pulchra_rebuild(self, model):
        """
        Calls pulchra for the first chain of the model, which is the one containing
        CA linkers

        :param model:   output list of [(PDBModel, modeled_domains, out_symseq)] 
                        as produced by Ranch, where for
                        each symmetric domain (only one if no symmetry) the
                        first chain is the one containing CA residues that must
                        be rebuilt
        :type model:    PDBModel
        
        :return m_reb:  model after the first chain was rebuilt by pulchra
        :type rebuilt:  PDBModel
        """

        m = model
        ch = m.takeChains([0])

        print('    Rebuilding with pulchra...')
        call = P.Pulchra(ch)
        ch_reb = call.run()

        m_reb = ch_reb.concat(m.takeChains(list(range(1,m.lenChains()))))

        m_reb.addChainId()
        m_reb['serial_number'] = N.arange(1,len(m_reb)+1) 

        return m_reb


    def restore_emb(self,emb_mod,emb_ch):
        '''
        Divides the emb_ch model into the original chains, as in emb_mod

        :param emb_mod: Model of the chains as they were initially embedded
        :type emb_mod:  PDBModel
        :param emb_ch:  Model of the chains to be separated in to the original ones
        :type emb_ch:   PDBModel
        '''

        l = 0   # Cummulative residue length of the chains
        full = B.PDBModel()
        
        # Go through each chain, obtain its length and take the corresponding
        # residues from emb_ch. Then renumber amino acids and concat to full
        for i in range(emb_mod.lenChains()):
            length = len(emb_mod.takeChains([i]).sequence())
            ch = emb_ch.takeResidues(list(range(l,l+length)))
            ch.renumberResidues()
            full = full.concat(ch)
            l += length

        return full


    ## WITHOUT GOING THROUGH ALL THE MODELS
    def extract_embedded(self,model,emb_mod,container_seq):
        """
        Extracts emb_mod from the first chain of the model, appends it at the end
        of the model and rebuilds first chain with pulchra

        :param model:   model to be processed
        :type model:    PDBModel
        :param emb_mod: PDBModel the embedded chain(s) that is located inside
                        of the first chain
        :type emb_mod:  PDBModel
        :param container_seq:   sequence of the specific domain with emb_mod inside
        :type container_seq:    str

        models = list of models to be processed
        emb_mod = PDBModel of the embedded CHAIN
        container_seq = sequence of the domain containing emb_mod
        symseq = sequence of the symmetric unit
        """
        processed = []

        m = model
       
        ch = m.takeChains([0])  # The first chain will contain the embedded chain(s)
        match = re.search(emb_mod.sequence(), ch.sequence())
        
        assert match, 'emb_mod is not in the first chain...'+'\n'+\
            emb_mod.sequence()+'\n'+ch.sequence()

        start, end = match.span()

        assert container_seq==ch.takeResidues(list(range(start-2, 
            start-2+len(container_seq)))).sequence(), 'emb_mod is not inside \
            container_seq...'

        emb_ch = ch.takeResidues(list(range(start,end)))
        atomstart = ch.resIndex()[start]
        atomend = ch.resIndex()[end]

        ch.remove(list(range(atomstart,atomend)))

        # Rebuilding the chain
        while ch.lenChains() > 1:
            ch.mergeChains(0)

        ch.renumberResidues()

        print('    Rebuilding with pulchra...')
        call_pulcura = P.Pulchra(ch)
        ch_reb = call_pulcura.run()

        # Concat ch_reb to the rest of the chains in m
        m_reb = ch_reb.concat(m.takeChains(list(range(1,m.lenChains()))))

        # Divide emb_ch into the original chains (after running ranch, the aa
        # are renumbered and the chains division is lost)
        emb_ch = self.restore_emb(emb_mod,emb_ch)

        # Concat the embedded chain(s) at the end
        m_reb = m_reb.concat(emb_ch)

        m_reb.addChainId()
        m_reb['serial_number'] = N.arange(1,len(m_reb)+1)

        return m_reb

    def embed_symmetric(self,j_doms, full_chains):
        """
        Embeds every symmetric unit (full_chain) from full_chains into every j_dom

        :param j_doms:  All the j_doms in which the symmetric chains will be
                        embedded
        :type j_doms:   list of PDBModels 
        """
        
        full_symmetric = B.PDBModel()
        emb_jsym = None
        ch = None
        # m = []
        
        for i in range(len(full_chains)):
            ch = full_chains[i]
            j_dom = j_doms[i]
            ch = self.extract_fixed(j_dom, ch)
            emb_jsym = R.embed(j_dom,ch)   # Symmetric unit embedded into j_dom

            assert emb_jsym.sequence() == j_dom.sequence()[:2] + ch.sequence() +\
                j_dom.sequence()[2:]

            full_symmetric = full_symmetric.concat(emb_jsym)

        emb_mod = ch
        symunit_ranch = emb_jsym.sequence()
        return full_symmetric, symunit_ranch, emb_mod

    def process_fullchain(self,chaini,model,out_symseq,bound_indexes):
        """
        Multi-purpose method that cleans the models obtained by ranch and saves
        some chain attributes
        
        It takes every symmetric unit by separate, and does the following:
        1.  Extracts the embedded (previously modeled) chains
        2.  Runs pulchra on the modeled chain
        3.  Adds symmetric unit (full_chain) and modeled domains dictionary to 
            chain properties
        4.  Replaces modeled domains in bound chains .domains attribute, and 
            removes said domains from the corresponding full_chain
        5.  Concatenates all of the full_chains into a full_model

        out_symseq = sequence of the entire symmetric unit produced by ranch, or
        the sequence of the full model if there is no symmetry. May contain
        embedded chains (emb_mod) inside certain domains (container_seq)

        :param chaini:  chain to be processed
        :type chaini:   multiprot.Chain object
        :param bound_indexes:   indexes of the chains bound to chaini
        :type bound_indexes:    list of int
        :param m:   model produced by Ranch (one of them)
        :type m:    PDBModel
        :param out_symseq:  sequence of the entire symmetric unit produced by
                            ranch, or the sequence of the full model if there is
                            no symmetry. It may contain embedded previously
                            modeled chains (emb_mod) inside certain domains
                            (container_seq)
        :type out_symseq:   str

        :return s:  number of symmetric units
        :type s:    int

        Other output/modification to external variables:
        This method appends every processed (emb_mod extracted and rebuilt with 
        pulchra) symmetric unit as a separate element to the self.full_chains 
        attribute, which is a list.

        """

        # if not full_model:
        # full_model = B.PDBModel()

        m = model[0]
        
        self.full_chains = []   # Reset variable to add new modeled chains

        # For every appearance of out_symseq sequence in the model
        s=0
        for match in re.finditer(out_symseq, m.sequence()):
            istart, iend = match.span()

            full_ch = m.takeResidues(list(range(istart, iend)))

            # If there is a previously modeled chain embedded somewhere
            if chaini.container_seq:
                # Extract the embedded chain(s) and rebuild with pulchra
                # embedded chain(s) end up at the end of the model
                full_ch = self.extract_embedded(full_ch, chaini.emb_mod, 
                    chaini.container_seq)
            else:
                # Only rebuild with pulchra
                full_ch = self.pulchra_rebuild(full_ch)

            # Add symmetric unit and modeled_domains dict to chain properties
            self.full_chains.append(full_ch)
            chaini.modeled_domains.append(model[1][s])

            self.replace_modeled(chaini,bound_indexes,s)

            s += 1

        return s


    def concat_full(self):
        """
        Concats self.full_chains into a single model
        """
        final = B.PDBModel()

        for m in self.full_chains:
            final = final.concat(m)

        final.addChainId()
        final['serial_number'] = N.arange(1,len(final)+1)

        print('Done.')

        return final

    def replace_jdoms(self,chainj):
        '''
        Replaces the models with new coordinates from chainj.new_domains 
        into chainj.domains, and deletes them from self.full_chains[0]

        This method is only for non-symmetric structures
        '''

        for i in range(len(chainj.domains)):
            # Take the domains with new coordinates in chainj.new_domains,
            # and if there isn't any take the original one
            j_dom = chainj.new_domains[i]
            if j_dom:
                chainj.domains[i] = j_dom
                self.full_chains[0] = self.extract_fixed(j_dom, self.full_chains[0])


    def create_full(self, i=0):
        """
        Method that will build the models through ranch and pulchra and return a
        list with 10 PDBModels
            MODIFY DESCRIPTION
        """

        chaini = self.CHAINS[i]

        # Take only 'num' number of models
        print('Chain %d' % (i+1))
        print('    Modeling with ranch...')
        models = self.call_ranch(chaini)[:self.num]

        model = models[0]   # take first PDBModel
        out_symseq = models[0][2]   # symmetric unit sequence... if there is no
                                    # symmetry, this will be the seq of the
                                    # entire model

        # Write each modeled chain to disk (optional!)
        # model[0].writePdb(os.path.join(self.dest,'chain%d.pdb' % i))
        
        # Find paired chains
        chaini.paired_to = self.find_paired(i)

        if not chaini.paired_to and len(self.CHAINS)>1:
            raise InputError('At least one of the chains used as input is not \
                bound to another one. If you want to model an individual chain \
                make a call to multiprot.py with only that chain')

        # Find indexes of bound chains
        bound_indexes = [key for key,value in chaini.paired_to.items()]
  
        s = self.process_fullchain(chaini,model,out_symseq, bound_indexes)

        chaini.modeled=True
        
        for j in bound_indexes:
            chainj = self.CHAINS[j]
            if not chainj.modeled:

                chainj.paired_to = self.find_paired(j)
                bound_j = [key for key, value in chainj.paired_to.items()]
                for k in bound_j:
                    if self.CHAINS[k].modeled:
                        # chaink = self.CHAINS[k]
                        # Add entire model to the first domain paired to j
                        pair = chainj.paired_to[k][0]
                        j_ind = chainj.names.index(pair[0])
                        # k_ind = chaini.names.index(pair[1])

                        if s>1: 
                            # There is symmetry, embed all symmetric units into
                            # domjs ... assume that chaink is chaini
                            j_doms = chainj.jdomains[j_ind]
                            emb_sym = self.embed_symmetric(j_doms, self.full_chains)
                            full_sym = emb_sym[0]
                            chainj.domains[j_ind] = full_sym
                            chainj.args["symtemplate"] = full_sym
                            chainj.args["symunit"] = emb_sym[1]
                            chainj.container_seq = emb_sym[1]
                            chainj.emb_mod = emb_sym[2]

                            # Delete chainj.args["fixed"] contents
                            chainj.args["fixed"] = []
                        else:
                            # Take domains from chainj.new_domains, and remove
                            # them from self.full_chains[0]
                            self.replace_jdoms(chainj)
                            
                            # Concat the full_chain to the first j_dom bound
                            # to chain k
                            jdom_new = chainj.domains[j_ind].concat(
                                self.full_chains[0])
                            chainj.args["fixed"].remove(chainj.domains[j_ind])
                            chainj.domains[j_ind] = jdom_new
                            chainj.args["fixed"].append(jdom_new)
                            # NOTE: Add jdom_new to chainj.args['chains'] dict??
                            # Not necessary so far since jdom_new is always the
                            # first chain, and that's the one taken if there is
                            # no chain specified in args['chains'] dict

                        self.create_full(j)
                        break

        return None

    def run(self):
        '''
        Calls methods to create chains and concatenate them
        '''
        self.create_full()
        return self.concat_full()


    def write_pdbs(self, models, dest, pref='mp'):
        '''
        Writes the pdbmodels to the specified destination
        '''

        f_out = [os.path.join(dest,pref+'_%02d.pdb' % i) for i in range(1,len(models)+1)]

        for i in range(len(models)):
            models[i].writePdb(f_out[i])

        return None    


#############
##  TESTING        
#############
import multiprot.parseChains as C
import multiprot.testing as testing

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
        self.builder1 = Builder(CHAINS1,args1.debug,args1.number,
            args1.destination)

        args4 = C.parsing(self.argstring4.split())
        CHAINS4 = C.create_chains(args4)
        self.builder4 = Builder(CHAINS4,args4.debug,args4.number,
            args4.destination)

        args5 = C.parsing(self.argstring5.split())
        CHAINS5 = C.create_chains(args5)
        self.builder5 = Builder(CHAINS5,args5.debug,args5.number,
            args5.destination)

        args2ch = C.parsing(self.argstring2ch.split())
        CHAINS2ch = C.create_chains(args2ch)
        self.builder2ch = Builder(CHAINS2ch,args2ch.debug,args2ch.number,
            args2ch.destination)

        args2chfixed = C.parsing(self.argstring2chfixed.split())
        CHAINS2chfixed = C.create_chains(args2chfixed)
        self.builder2chfixed = Builder(CHAINS2chfixed,args2chfixed.debug,
            args2chfixed.number,args2chfixed.destination)

        args3ch = C.parsing(self.argstring3ch.split())
        CHAINS3ch = C.create_chains(args3ch)
        self.builder3ch = Builder(CHAINS3ch,args3ch.debug,args3ch.number,
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
