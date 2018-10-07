"""
Script that contains the higher level implementation of the ranch wrapper,
to handle single and multiple chain scenarios

"""


##### FOR THE EARLY IMPLEMENTATION OF MULTIPROT, WE WILL ASSUME ONLY ONE CHAIN
##### IN THE INPUT

import ranch as R
import pulchra as P
import tempfile
import os, re, operator
import biskit as B
import biskit.tools as T
import numpy as N

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
        self.CHAINS = chains    # Original chains and PDBModels from input

        self.dest = dest
        self.debug = debug

    def find_paired(self, i):
        """
        Finds to which chain(s) is chain i paired to
        
        :param i: Index of the chain to be compared against
        :type i: int

        :return paired_to:  dictionary with j:[[(namei),(namej)],[(namei),(namej)]] 
                            elements where j is the index of the chain paired to
                            chain 'i', namei is the name of the domain in chain i, 
                            and namej the name of the paired domain in chain j
        """

        paired_to = []

        chain_idxs = list(range(len(self.CHAINS)))
        chain_idxs.remove(i)

        for j in chain_idxs:
            for key, value in self.CHAINS[j].args['chains']:
                for pdb, chainid in self.CHAINS[i].args['chains']:
                    if key==pdb and value!=chainid:
                        if j in paired_to:
                            paired_to[j].append([(pdb,chainid),(key,value)])
                        else:
                            paired_to[j] = [[(pdb,chainid),(key,value)]]

                    # If same pdb and chain are referenced in different chains
                    # in args['chains'] parameter, they are not bound

        return paired_to

    def pulchra_rebuild(self, models):
        """
        Calls pulchra for the chains with CA residues, and cleans the rebuilt
        models

        :param models:  output list of [(PDBModel, modeled_domains, out_symseq)] 
                        as produced by Ranch, where for
                        each symmetric domain (only one if no symmetry) the
                        first chain is the one containing CA residues that must
                        be rebuilt
        :type models:   list
        :return rebuit: list of 10 PDBModels after the chains with CA residues
                        were rebuilt by pulchra
        :type rebuilt:  list
        """
        rebuilt = []

        # create temporary folder to write models
        self.tempdir = tempfile.mkdtemp('', self.__class__.__name__.lower() + '_')
        # create names to write pdbs into
        pdb_paths = [os.path.join(self.tempdir, 'm%02d' % i) for i in range(1,11)]

        for i in range(len(models)):

            # Rebuild only the chains with CA, which is the first one plus its
            # repetitions in symmetric units

            full = B.PDBModel()
            m = models[i]
            to_reb_seq = m.takeChains([0]).sequence()
            to_reb_len = len(m.takeChains([0]))

            for j in range(m.lenChains()):
                chain = m.takeChains([j])

                if chain.sequence() == to_reb_seq and len(chain) == to_reb_len:
                    chname = pdb_paths[i]+'_ch'+str(j)+'.pdb'
                    chain.writePdb(chname)

                    call = P.Pulchra(chname)
                    rebuild = call.run()
                    full = full.concat(B.PDBModel(chname[:-3]+'rebuilt.pdb'))

                else:
                    full = full.concat(chain)
        
            full.addChainId()
            full['serial_number'] = N.arange(1,len(full)+1)
            rebuilt.append(full)

        return rebuilt

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
            if N.all(dom.xyz[lowdom:highdom+1] == full.xyz[lowfull:highfull+1]):
                atom_start = full.resIndex()[start]
                atom_end = full.res2atomIndices([end-1])[-1] + 1
                full.remove(list(range(atom_start,atom_end)))
                break

        return full

    def replace_modeled(self, chaini, bound_indexes, s):
        """
        MODIFY
        Replaces the coordinates of the domains already modeled in a previous
        chain in the 'domains' attribute of chain j
        
        Replaces only the domains that are bound to other chains.

        In chain j it leaves only the part of the domain belonging to that chain (j_dom)
        And removes j_dom from the (already modeled) chain i
        
        :return None:   The method will only replace the corresponding domains 
                        in self.CHAINS
        """
        
        modeled_domains = chaini.modeled_domains[s]
        # For every bound chain
        for j in bound_indexes:
            chainj = self.CHAINS[j]
            # For every pair of coupled domains [(namei),(namej)]
            for pair in chaini.paired_to[j]:
                # Find the index of the domain in chains i and j
                i_ind = chaini.names.index(pair[0])
                j_ind = chainj.names.index(pair[1])
                
                # Replace domain with new coordinates in chain j, and fix
                # Replace with only the part belonging to chain j
                j_id = pair[1][1]
                modeled_dom = modeled_domains[i_ind]
                j_chainind = chainIndex(modeled_dom,j_id)
                j_dom = modeled_dom.takeChains(j_chainind)
                
                # Only replace in 'domains' atribute the j_dom from the first
                # symmetric unit
                if s=0:
                    self.CHAINS[j].domains[j_ind] = j_dom
                    self.CHAINS[j].args["fixed"].append(j_dom)

                # Remove j_dom from chain i
                chaini.full_chain[s] = R.extract_fixed(j_dom, chaini.full_chain[s])
            
                ###########
                ## CHECK CHAINS ONE BY ONE, EXTRACT THOSE THAT ARE EQUAL TO J_DOM

        return None

    def chainIndex(self,model,chain_id):
        """
        Take as input the chain id of model and return its index
        """
        mask_chain = model.maskFrom('chain_id', chain_id)
        i_mask_chain = N.nonzero(mask_chain)[0]
        return model.atom2chainIndices(i_mask_chain)[0]

    def repalce_symmetric(self,chaini,i_ind):
        """
        Embeds every symmetric unit into every domain i
        """
        full_model = B.PDBModel()
        emb_sym = None
        m = []
        
        for i in range(len(chaini.full_chain)):
            ch = chaini.full_chain[i]
            m.append(chaini.modeled_domains[i][i_ind])
            emb_sym = R.embed(m[i],ch)   # Symmetric unit embedded into m
            full_model = full_model.concat(emb_sym[i])

        symunit = emb_sym.sequence()
        return full_model, m, symunit

    def cleanup(self):
        """
        Delete temporary files
        """
        if not self.debug:
            T.tryRemove(self.tempdir, tree=True)

    def ranch_pulchra(self,chaini):
        """
        Builds model with ranch, and rebuilds CA chains with pulchra
        """

        # Model with ranch
        call = R.Ranch(*chaini.domains, **chaini.args, debug=self.debug)
        models = call.run()

        # Rebuild CA chains with pulchra
        pdbmodels = [p[0] for p in models]
        rebuilt = pulchra_rebuild(pdbmodels)

        # Replace rebuilt models
        for i in range(len(models)):
            models[i][0] = rebuilt[i]

        return models

    def process_symunits(self,chaini, out_symseq, bound_indexes):


    def model(self, i=0, full_model = None, emb=[]):
        """
        Method that will build the models through ranch and pulchra and return a
        list with 10 PDBModels
            MODIFY DESCRIPTION
        """
        if not full_model:
            full_model = B.PDBModel()

        chaini = self.CHAINS[i]

        models = ranch_pulchra(chaini)

        m = models[0][0]   # take first PDBModel
        out_symseq = models[0][2]   # symmetric unit sequence... if there is no
                                    # symmetry, this will be the seq of the
                                    # entire model
        
        # Find paired chains
        chaini.paired_to = find_paired(i)

        # Find indexes of bound chains
        bound_indexes = [key for key,value in chain.paired_to.items()]

        # For every appearance of out_symseq sequence in the model
        s=0
        for match in re.finditer(out_symseq, m.sequence()):
            istart, iend = match.span()

            # Add symmetric unit and modeled_domains dict to chain properties
            chaini.full_chain.append(m.takeResidues(list(range(istart, iend))))
            chaini.modeled_domains.append(models[0][1][s])

            # Remove chain ids, so there are no conflicts with ids from chains 
            # to be modeled
            chaini.full_chain[s].atoms['chain_id'] = ['#'] * len(chaini.full_chain)
            
            replace_modeled(chaini,bound_indexes,s)
            full_model = full_model.concat(chaini.full_chain[s])
            s += 1

        chaini.modeled=True

        for j in bound_indexes:
            chainj = self.CHAINS[j]
            if not chainj.modeled:
                chainj.paired_to = find_paired(j)
                bound_j = [key for key, value in chainj.paired_to.items()]
                for k in bound_j:
                    if self.CHAINS[k].modeled:
                        chaink = self.CHAINS[k]

                        # Add entire model to the first domain paired to j
                        pair = chainj.paired_to[k][0]
                        j_ind = chainj.names.index(pair[0])
                        k_ind = chaini.names.index(pair[1])

                        if s>1: 
                            # There is symmetry, embed all symmetric units into
                            # domjs ... assume that chaink is chaini
                            rep_sym = replace_symmetric(chaink,k_ind)
                            full_model = rep_sym[0]
                            chainj.domains[j_ind] = full_model
                            chainj.args["symtemplate"] = full_model
                            chainj.args["symunit"] = rep_sym[2]
                            emb.append(rep_sym[1])
                        else:
                            j_dom = chainj.domains[j_ind].concat(full_model)
                            self.CHAINS[j].domains[j_ind] = j_dom

                model(j,full_model,emb)

        return full_model, emb

    def extract_symmetric(full, emb):
        """
        Extracts each emb domain from full_model, backwards
        """
        inds = sorted(list(range(len(emb))), reverse=True)

        for i in inds:
            for m in emb[i]:
                # Remove the embedded sequences and add them at the end of full
                full = Ranch.extract_fixed(m,full)
                full = full.concat(m)

        return full

    def build(self):
        result = model()
        full_model = result[0]
        emb = result[1]

        if emb:
            full_model = extract_symmetric(full_model, emb)

        full_model.renumberResidues()
        full_model.addChainId()
        full_model['serial_number'] = N.arange(1,len(full_model)+1)

        return full_model

#############
##  TESTING        
#############
import testing
import multiprot as mp

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