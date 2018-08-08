"""
Script that contains the higher level implementation of the ranch wrapper,
to handle single and multiple chain scenarios

"""


##### FOR THE EARLY IMPLEMENTATION OF MULTIPROT, WE WILL ASSUME ONLY ONE CHAIN
##### IN THE INPUT

import ranch as r
import pulchra as p
import tempfile
import os

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

    def build(self):
        """
        Method that will build the models and return a return a list with 10
        PDBModels
        """

        # Call ranch for every chain
        for i in range(len(chains)):
            if chains[i][5] == False:   # If it is still not modeled
                call = r.Ranch(*chains[i][0], chains=chains[i][1], 
                    symmetry=args.symmetry, fixed = chains[i][2], 
                    symtemplate = chains[i][3], pool_sym=chains[i][4],
                    debug=self.debug)
                models = call.run()
                chains[i][4] = True
                
        ## ***** code for multiple chains goes here

        # create temporary folder to write models
        self.tempdir = tempfile.mkdtemp('', self.__class__name__.lower() + '_')
        # create names to write pdbs into
        pdb_paths = [os.path.join(tdir, 'm%02d.pdb' % i) for i in range(1,11)]

        # Write pdb, call pulchra and move rebuilt pdb to destination
        for i in range(len(models)):
            models[i].writePdb(pdb_paths[i])
            call = p.Pulchra(pdb_paths[i])
            rebuild = call.run()
            os.rename(pdb_paths[i][:-3]+'rebuilt.pdb', 
                os.path.join(self.dest, 'm%02d.pdb' % i))

    def cleanup(self):
        """
        Delete temporary files
        """
        if not self.debug:
            T.tryRemove(self.tempdir, tree=True)