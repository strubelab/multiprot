"""
Script that contains the higher level implementation of the ranch wrapper

"""

#RG: I think this should also be a class with methods, defaults and all
#RG: the actual script should then be very small and only deal with command line parsing and argument cleanup

import argparse
import biskit as b
import random
import ranch as r

# If type=divide does not work, try action='append_const'
#RG: no big deal but instead of `string`, which looks like a type definition, convention is `s`
def divide(string):  
   return tuple(string.split(':'))
   # if the string is not properly formatted
   # raise argparse.ArgumentTypeError(msg)

# Supported symmetries for --symmetry argument
supp_sym = ['p'+str(i) for i in range(1,20)]
supp_sym += ['p'+str(i) for i in range(22,132,10)] + ['p222']

# Add positional and optional arguments
# NOTES: try formatter_class=argparse.MetavarTypeHelpFormatter
#        if no argument is given, show help (and possibly raise error)
#        look at customizing file parsing
#        look at exiting methods

parsero = argparse.ArgumentParser(usage='', 
   description='''Build multiple chain-proteins. The arguments can also
   be read from a file, in which case the file name must have the @ prefix''')

parsero.add_argument('--chain', '-c', action='append', nargs='+', type=divide, 
   help='Add a new chain to the model', required=True)

parsero.add_argument('--split', '-spl', default=None, 
   help='Split one or more chains from a given PDB')

parsero.add_argument('--symmetry', '-sym', default='p1',
   help='What kind of symmetry do you wish to have in your molecule. Supported\n\
   symmetries are: p1, p2, …, p19 (nineteen-fold), p22, p32, p42, p52, p62,\n\
   …, p122, p222.', choices=supp_sym)

parsero.add_argument('--symtemplate', '-t', default=[], action='append' 
   help='Which domain will be the symmetry core, in case of symmetry other\
   than p1 specified')

parsero.add_argument('--overall_sym', '-o', default='mix', 
   choices=['mixed', 'm', 'symmetry', 's', 'asymmetry', 'a'], 
   help='Specify the overall symmetry of the molecules to be produced, i.e. \
   all symmetric [s], all asymmetric [a] or mixed. [m]')

parsero.add_argument('--fixed', '-f', default=[], nargs='*'
   help='Specify one or more domains to be fixed in their original coordinates.')

parsero.add_argument('args', nargs=argparse.REMAINDER)

#argument_default = argparse.SUPPRESS

args = parsero.parse_args()

# parsero.Namespace() ?

# vars(args)  returns dictionary with attributes

## Create PDBModels

chains = []

for i in range(len(args.chain)):    # For each chain
   rdomains = []
   rchains = {}
   rfixed = []
   rsymtemp = None

   for j in range(len(args.chain[i])):    # For each component of the chain
      if args.chain[i][j][0][-4:]=='.pdb':
         pdb = b.PDBModel(args.chain[i][j][0])
            if len(args.chain[i][j])==2:
               rchains[pdb] = args.chain[i][j][1]
            if args.chain[i][j][0] in args.fixed:
               rfixed.append(pdb)
            if args.chain[i][j][0] in args.symtemplate:
               rsymtemp = pdb    # THERE CAN ONLY BE ONE

         rdomains.append(pdb)
         continue
      rdomains.append(args.chain[i][j][0])

   # chains[i] = (domains/linkers, chains dict, symtemplate, already modeled)
   chains.append((rdomains, rchains, rfixed, rsymtemp, False))

random.shuffle(chains)

for chain in chains:
   if chain[2] == False:
      call = r.Ranch(*chain[0], chains=chain[1], symmetry=args.symmetry, 
         symtemplate=args.symtemplate, overall_sym=args.overall_sym,)



args=parsero.parse_args('--chain dom1.pdb:A ECSEIVIERECSEIVIER dom2.pdb --chain dom1.pdb:B ECSEIVIERECSEIVIER dom4.pdb'.split())
