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





############## RUN NEW_RANCH_5 ##################

## EXAMPLE 1

import multiprot.new_ranch_5 as r
import biskit as b

dom1 = b.PDBModel('/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/2z6o_mod.pdb')
dom2 = b.PDBModel('/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/Histone_H3.pdb')
call = r.Ranch(dom1,'GGGGGGGGGG',dom2)
models = call.run()

## EXAMPLE 4

import multiprot.new_ranch_5 as r
import biskit as b

domAB1 = b.PDBModel( "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
domAB2 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', domAB2, chains = {domAB1:'A', domAB2: 'B'})
models = call.run()

## EXAMPLE X

import multiprot.new_ranch_5 as r
import biskit as b

domAB1 = b.PDBModel( "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
domAB2 = domAB1.clone()
domAB3 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', domAB2, 'GGGGGGGGGGGGGGGGGGGG', domAB3, chains = {domAB1:'A', domAB2: 'A', domAB3:'B'})
models = call.run()

## EXAMPLE 5
## HAVE TO FIX CLEANING OF SEQUENCE POST-RUN

import multiprot.new_ranch_5 as r
import biskit as b

domAB1 = b.PDBModel( "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
domAB2 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', domAB2, chains = {domAB1:'A', domAB2: 'B'}, 
    symmetry='p2', symtemplate=domAB2)
models = call.run()

## Extra

import multiprot.new_ranch_5 as r
import biskit as b

domAB1 = b.PDBModel( "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
dom2 = b.PDBModel('/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/Histone_H3.pdb')
domAB2 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', dom2, 'GGGGGGGGGGGGGGGGGGGG', domAB2, 
    chains = {domAB1:'A', domAB2: 'B'})
models = call.run()


## Save PDBModels
filenames = ['/Users/guzmanfj/Documents/Stefan/multiprot/multiprot/tests/m' + str(i) + '.pdb' for i in range(10)]

for i in range(len(models)):
    models[i].writePdb(filenames[i])
