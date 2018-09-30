# multiprot
Model proteins with multiple domains

This project wraps and adds functionality to two different programs:
- Ranch from the [Ensemble Optimization Method](https://www.embl-hamburg.de/biosaxs/eom.html)
- [Pulchra](http://www.pirx.com/pulchra/index.shtml)

You can produce models combining structured domains, providing the corresponding atomic coordinates in PDB format, and random coil or "linker" domains.

This is an early version which only offers support for models where the linker segments are located on a single chain.


INSTRUCTIONS FOR THE PACKAGE TO USE IN THE COMMAND LINE

Go to the package directory and run the `multiprot.py` module with the required arguments.
To see a full list of arguments run `python multiprot.py -h`. Here you can find a few examples, where all the file names are relative to the `multiprot/multiprot` package folder:

## Example 1
Build simple single-chain protein joining two structured domains (2z6o.pdb and histone.pdb) with a disordered linker (GGGGGGGGGG), and saving the resulting models to example1/ directory
```
$ python multiprot.py --chain testdata/2z6o.pdb GGGGGGGGGG testdata/histone.pdb --destination example1
```

## Example 2
Build model joining domAB1.pdb (chain A) with domAB1.pdb (chain B) through a Gly linker and saving the resulting models to example2/
```
$ python multiprot.py --chain testdata/domAB1.pdb:A GGGGGGGGGG testdata/domAB1.pdb:B --destination example2
```

## Example 3
Symmetric protein taking a dimer as a symmetry core, and adding a linker and another dimeric domain (four chains total)
```
python multiprot.py --chain testdata/domAB1.pdb GGGGGGGGGG testdata/domAB1.pdb:A --symmetry p2 --symtemplate testdata/domAB1.pdb --poolsym s --destination example3
```



((will continue...))

freezer - include Python

pip install -r requirements.txt

pip install ./downloads/SomeProject-1.0.4.tar.gz