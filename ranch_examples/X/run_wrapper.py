## MADE WITH RANCK WRAPPER

import multiprot.new_ranch_5 as r
import biskit as b

domAB1 = b.PDBModel( "/Users/guzmanfj/Documents/Stefan/multiprot/ranch_examples/4/dom1_AB.pdb")
domAB2 = domAB1.clone()
domAB3 = domAB1.clone()
call = r.Ranch(domAB1, 'GGGGGGGGGGGGGGGGGGGG', domAB2, 'GGGGGGGGGGGGGGGGGGGG', domAB3, 
    chains = {domAB1:'A', domAB2: 'A', domAB3:'B'})
models = call.run()     # List with 10 PDBModels


## Save PDBModels
filenames = ['/Users/guzmanfj/Documents/Stefan/multiprot/multiprot/tests/m' + str(i) \
    + '.pdb' for i in range(10)]

for i in range(len(models)):
    models[i].writePdb(filenames[i])
