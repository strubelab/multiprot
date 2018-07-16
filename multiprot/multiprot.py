"""
Script that contains the higher level implementation of the ranch wrapper

"""

import argparse

# If type=divide does not work, try action='append_const'
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
   help='What kind of symmetry do you wish to have in your molecule. Supported\
   symmetries are: p1, p2, …, p19 (nineteen-fold), p22, p32, p42, p52, p62,\
   …, p122, p222.', choices=supp_sym)

parsero.add_argument('--symtemplate', '-t', default=None, 
   help='Which domain will be the symmetry core, in case of symmetry other\
   than p1 specified')

parsero.add_argument('--overall_sym', '-o', default='mix', 
   choices=['mixed', 'm', 'symmetry', 's', 'asymmetry', 'a'], 
   help='Specify the overall symmetry of the molecules to be produced, i.e. \
   all symmetric [s], all asymmetric [a] or mixed. [m]')

parsero.add_argument('--fixed', '-f', default=None, 
   help='Specify a domain to be fixed in its original coordinates.')

parsero.add_argument('args', nargs=argparse.REMAINDER)

#argument_default = argparse.SUPPRESS

args = parsero.parse_args()

# parsero.Namespace() ?

vars(args)  # returns dictionary with attributes

