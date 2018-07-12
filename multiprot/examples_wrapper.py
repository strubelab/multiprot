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





############## RUN NEW_RANCH_6 ##################

## EXAMPLE 1

import multiprot.new_ranch_6 as r
import biskit as b

dom1 = b.PDBModel(
    '/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/2z6o_mod.pdb')
dom2 = b.PDBModel(
    '/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/Histone_H3.pdb')
call = r.Ranch(dom1,'GGGGGGGGGG',dom2)
models = call.run()

## EXAMPLE 4

import multiprot.new_ranch_6 as r
import biskit as b

domAB1 = b.PDBModel(
    "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
domAB2 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', domAB2, 
    chains = {domAB1:'A', domAB2: 'B'})
models = call.run()

## EXAMPLE 5

import multiprot.new_ranch_6 as r
import biskit as b
 
domAB1 = b.PDBModel(
    "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
domAB2 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', domAB2, chains = {domAB2: 'A'}, 
    symmetry='p2', symtemplate=domAB1, overall_sym='symmetry')
models = call.run()

## EXAMPLE 7

import multiprot.new_ranch_6 as r
import biskit as b

# No need to specify chains dict entry for multichain domain if it is the symtemplate
dom2 = b.PDBModel(
    '/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/Histone_H3.pdb')
domAB = b.PDBModel(
    "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
linker = 'GGGGGGGGGGGGGGGGGGGG'
call = r.Ranch(dom2, linker, domAB, linker, dom2, 
    symmetry='p2', symtemplate=domAB, overall_sym='mixed')
models = call.run()

## EXAMPLE 10

import multiprot.new_ranch_6 as r
import biskit as b

# If chain to be taken is not specified it defaults to the first chain

# It is not necessary to clone a multichain domain and specify the chain to be taken
# if the chain will be the same for all the appearances of the domain

domAB1 = b.PDBModel(
    "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
domAB2 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', domAB2, 'GGGGGGGGGGGGGGGGGGGG', domAB2, 
    chains = {domAB2:'B'})
models = call.run()

## Extra

import multiprot.new_ranch_6 as r
import biskit as b

domAB1 = b.PDBModel(
    "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
dom2 = b.PDBModel(
    '/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/1/Histone_H3.pdb')
domAB2 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', dom2, 'GGGGGGGGGGGGGGGGGGGG', domAB2, 
    chains = {domAB1:'A', domAB2: 'B'})
models = call.run()


## Save PDBModels
filenames = ['/Users/guzmanfj/Documents/Stefan/multiprot/multiprot/tests/m' + str(i) + '.pdb' for i in range(10)]

for i in range(len(models)):
    models[i].writePdb(filenames[i])

