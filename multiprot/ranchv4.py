"""
Author: Francisco Javier Guzman

Last edited: 20/feb/18

"""

"""

A ranch wrapper to generate a pool of n independent models based upon a
sequence and structural information. For multi-domain proteins where high-
resolution structures for individual domains are available, such files (e.g. 
PDB) can be used as rigid-body domains/subunits and/or as an aid to define 
distance constraints during model generation. For proteins expected to be 
intrinsically unfolded, no rigid bodies are used and random configurations 
of the alpha-carbon trace are created based upon the sequence. Crystallographic 
symmetry can anso be applied (P1, P2, Pn) and requires a high-resolution
multichain/oligomerized PDB file as input, or through specification of a
potential oligomerization interface via a set of user defined distance
constraints.

Application:
Reference:


INPUT EXAMPLES FOR RANCH

Parameters:

sequence.seq 	file with the sequence of the model to be built. The order of the
				units and linkers is important
-s <VALUE>		symmetry for the molecule

-x dom.pdb		files with the domain of known structure. It can have one or
				multiple chains (units)
-f <yes|no>		maintain the original coordinate position of this rigid body for 
				each model generated
-o <yes|no>		input PDB file contains a multichain interface

-f <SUFFIX>		Suffix name for the each model/PDB file created

-w <PATH>		Path to the folder where the models are to be saved

-c <chainlgth>	Length for each chain in the final model

	
	Different rigid bodies (RBs) are labeled with letters A, B, C etc., \
	and individual chains belonging to the same RB are labeled A1, A2, etc.

	Example1: 2domains_exm_p2
	A1~B x 2
	- 2 rigid bodies/domains (A, B)
	- 1 fixed rigid body (RB) with 2 identical chains (A1A2), and 1 \
		'loose' RB (B)
	- Total: 2 identical chains with domain A1-linker-domain B

	ranch sequence.seq -s=p2 -q=10 -x=dom1.pdb -x=dom2.pdb -f=yes -f=no 
	--filesuffix=2d_exmP2 -o=yes -o=no -w=2domains_exm_p2 -i
	
	sequence.seq:
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	SDLIEGRGIP
	MISLIAALAVDRVIGMENAMPWNLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRKNIILSSQPGTDDRVTWVKSVDE
	AIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFEILERR

	Example2: 2domAB_seqmod_p1
	A1A2~B
	- 2 RBs (A, B)
	- 1 fixed RB with 2 identical chains (A1A2), and 1 'loose' RB (B)
	- Total: 2 chains, A1 and A2-linker-B

	$ ranch sequence_mod.seq -q=10 -x=dom1_AB.pdb -x=dom2.pdb -f=yes -f=no 
	--filesuffix=AB_seqmod_p1 -o=yes -o=no -w=2domAB_seqmod_p1 -i
	
	sequence_mod.seq:
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	SDLIEGRGIP
	MISLIAALAVDRVIGMENAMPWNLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRKNIILSSQPGTDDRVTWVKSVDE
	AIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFEILERR
	
	Example3: 2dom_2chains
	A1~B1B2~A2
	- 2 identical RBs/domains (A, B)
	- 1 fixed RB with 2 identical chains (A1A2), and 1 \
		'loose' RB with two identical chains (B1B2)
	- Total: 2 chains, A1-linker-B1, B2-linker-A2
	
	$ ranch sequence_mod4.seq -q=10 -x=dom1_A.pdb -x=dom1_AB.pdb -x=dom1_B.pdb -f=yes 
	-f=no -f=yes --filesuffix=2d_2ch -o=no -o=yes -o=no -w=2dom_2chains -i
	
	sequence_mod4.seq:
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	SDLIEGRGIPSDLIEGRGIPSDLIEGRGIPSDLIEGRGIP
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	SDLIEGRGIPSDLIEGRGIPSDLIEGRGIPSDLIEGRGIP
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK

	Example4: 2dom_linear
	A1~B~A2
	- 2 RBs (A,B)
	- 1 fixed RB with 2 identical chains (A1A2), and 1 'loose' RB (B)
	- Total: 1 chain, A1-linker-B-linker-A2
	
	$ ranch sequence_mod2.seq -q=10 -x=dom1_A.pdb -x=dom2.pdb -x=dom1_B.pdb -f=yes 
	-f=no -f=yes --filesuffix=2dom_linear -o=no -o=no -o=no -w=2dom_linear -i

	sequence_mod2.seq
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	SDLIEGRGIPSDLIEGRGIP
	MISLIAALAVDRVIGMENAMPWNLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRKNIILSSQPGTDDRVTWVKSVDE
	AIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFEILERR
	SDLIEGRGIPSDLIEGRGIPSDLIEGRGIP
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK

	3dom_linear
	A1~B~A1~C
	- 3 RBs (A,B,C)
	- 1 fixed RB with 2 identical chains (A1A2), and 2 'loose' RB (B,C)
	- Total: 1 chain, A1-linker-B-linker-A2-linker-C
	
	$ ranch sequence_mod3.seq -q=10 -x=dom1_A.pdb -x=dom2.pdb -x=dom1_B.pdb 
	-x=histone.pdb -f=yes -f=no -f=yes -f=no --filesuffix=3dom_linear -o=no -o=no 
	-o=no -o=no -w=3dom_linear -i

	sequence_mod3.seq
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	SDLIEGRGIPSDLIEGRGIP
	MISLIAALAVDRVIGMENAMPWNLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRKNIILSSQPGTDDRVTWVKSVDE
	AIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFEILERR
	SDLIEGRGIPSDLIEGRGIPSDLIEGRGIP
	MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHN
	MLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALD
	VVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK
	MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK
	TDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA


"""

import biskit as B
import numpy as N
import tempfile, os

from biskit.exe.executor import Executor
import biskit.tools as t

class Ranch( Executor ):
	"""
	Ranch wrapper

	Usage
	=====

	NOTES:
	## Input (sequence and proteins) as python objects or file names?
	## List of domains, fixed and multichain values for each rigid body
	## Return list with generated PDBModels
	## Add more parameters for full functionality of ranch
	## Check for '-' characters

	>>> prot = Ranch( sequence, sym, doms, fixed, multich, filesuff, savepath)

	Ranch takes none or several PDBModel instances as input for rigid bodies and 
	the sequence of the final model and returns a list containing 10 models with
	the integrated structures in the order of the given sequence.

	... more ranch description (link to manual?)...

	"""

	def __init__(self, sequence=None, s=None, x=None, f=None, o=None,
		filesuff=None, **kwargs):
		"""
		:param sequence: sequence of the final protein model, with rigid \
				bodies and linkers
		:type sequence: str
		
		:param kw: additional key=value parameters for ranch (See ranch \
			manual). Same parameters as in the normal ranch usage (with \
			the exception of 'filesuff') but without the '-' character \
			and instead of providing the same parameter multiple times \
			for every domain, input the argument once with a tuple containing \
			the values for each domain.
			::
			s=<VALUE>		symmetry for the molecule ('p1', 'p2', etc.)

			x=(dom1,dom2,)	PDBModels of the domains/subunits that you \
							wish to define as rigid bodies.

			f=('yes','no',) maintain the original coordinate position of \
			default='no'	this rigid body for each model generated.
			
			o=('yes','no',) select 'yes' if the input PDBModel contains a \
			default='no'	multichain interface (symmetric or asymmetric) \
							that can be used as the structural core of each \
							model. Select 'no' if the PDBModel defines a \
							monomeric subunit/domain

			filesuff=<SUFFIX>	Suffix name for the each model/PDB file created ... necessary?
			default='eom'

			w <PATH>		Path to the folder where the models are to be saved
			default='.'

		"""

		# Create temporary folder for pdbs and sequence
		tempdir = tempfile.mktemp( '', self.__class__.__name__.lower() + '_', t.tempDir() )

		# Create temporary folder for models
		self.tempmodels = tempfile.mktemp( '', 'models_', tempdir )

		self.sequence = sequence
		self.f_seq = tempfile.mktemp('_sequence.seq', '', tempdir)

		self.s = s
		self.x = x
		self.f = f
		self.o = o
		self.filesuff = filesuff
		
		# Generate by default 10 models, with no intensities
		args = self.f_seq + ' -q=10 -i'

		if self.s:
			args = args + ' -s=%s' % self.s

		if self.x:
			self.pdbs_in = ()
			for dom in self.x:
				# Best way to name? Maybe it doesnt matter
				self.pdbs_in = self.pdbs_in + (tempfile.mktemp(
					dom.sourceFile()[-8:], '', tempdir),)

			args = args + ' -x=%s' * len(self.pdbs_in) % self.pdbs_in

		if self.f:
			args = args + ' -f=%s' * len(self.f) % self.f

		if self.o:
			args = args + ' -o=%s' * len(self.o) % self.o

		if self.filesuff:
			args = args + ' --filesuffix=%s' % (self.filesuff,)

		args = args + ' -w=%s' % (self.tempmodels)

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
		