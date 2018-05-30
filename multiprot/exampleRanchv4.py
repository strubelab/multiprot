## TEST RANCHV3

import biskit as b
import ranchv3 as r

dom1 = b.PDBModel("ranch_example/2z6o_mod.pdb")
dom2 = b.PDBModel("ranch_example/Histone_H3.pdb")

with open("ranch_example/full.fasta") as f:
	seq=f.readlines()

del(seq[0])

seq = ''.join(seq).replace('\n','')

call = r.Ranch(sequence=seq, s='p1', x=(dom1,dom2), f=('yes','no'), o=('no','no'), filesuff='test', w='ranch_example/models')

## TEST RANCHV4
## cwd = Documents/Stefan/multiprot/multiprot

import biskit as b
import ranchv4 as r

dom1 = b.PDBModel("ranch_example/2z6o_mod.pdb")
dom2 = b.PDBModel("ranch_example/Histone_H3.pdb")

with open("ranch_example/full.fasta") as f:
	seq=f.readlines()

del(seq[0])

seq = ''.join(seq).replace('\n','')

call = r.Ranch(sequence=seq, s='p1', x=(dom1,dom2), f=('yes','no'), o=('no','no'), filesuff='test')

models = call.run()		# should return a list with 10 models
