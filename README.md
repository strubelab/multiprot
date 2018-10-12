# multiprot
Model proteins with multiple domains

This project wraps and adds functionality to two different programs:
- Ranch from the [Ensemble Optimization Method](https://www.embl-hamburg.de/biosaxs/eom.html)
- [Pulchra](http://www.pirx.com/pulchra/index.shtml)

You can produce models combining structured domains, providing the corresponding atomic coordinates in PDB format, and random coil or "linker" domains.

This is an early version which only offers support for models where the linker segments are located on a single chain.

## multiprot installation

### Install Ranch helper application

Ranch is part of the ATSAS software package from the EMBL Hamburg (BioSax) team. It was originally a stand-alone program which has now been merged into the [EOM package](https://www.embl-hamburg.de/biosaxs/eom.html). 

1. Create a user account at [ATSAS website account registration](https://www.embl-hamburg.de/biosaxs/atsas-online/register.php)
2. Download ATSAS package for your architecture from [ATSAS download area](https://www.embl-hamburg.de/biosaxs/atsas-online/download.php)
3. Install package. E.g. for Debian or Ubuntu:
    ```sh
    sudo apt install ./ATSAS-2.8.4-1_amd64.deb
    ```
    This will, among other tools, put the `ranch` command into your search path. Open a terminal and verify that the program is correctly installed. Typing the `ranch` command should open the interactive ranch prompt:
    ```sh
    ~> ranch
    *******  ------------------------------------------------------  *******
    *******     RANCH - version 2.2 - (r10552)                    ********
    ...
    ```

### Install Pulchra helper application

to do: describe

### Install multiprot Python package

For developers, installation directly from github should work like this: 

1. Download latest multiprot version (this will create a folder `multiprot` within the current directory):
   ```sh
   git clone https://github.com/StruBE-KAUST/multiprot.git
   ```
2. set up and start virtual environmnent (This step is optional. Skip if you prefer installing multiprot system-wide):
   ```sh
   cd multiprot
   python3 -m venv venv
   source venv/bin/activate
   ```
   At this point, it's a good idea to update your `pip` tool (only affects the virtual environment):
   ```sh
   pip install --upgrade pip
   ```
3. multiprot depends on the python 3 version of [biskit](https://github.com/graik/biskit) which is not yet available on pypi so we have to install it by hand:
   ```sh
   git clone https://github.com/graik/biskit.git biskit
   pip3 install -r biskit/requirements.txt
   pip3 install -e biskit
   ```
4. Ensure multiprot is in the virtual environment's python path. Easiest is to create a symbolic link because that avoids any issues with having two different copies of the package flying around on your system. From within the base `multiprot` folder (as before) run:
   ```sh
   ln -s multiprot venv/lib/venv/lib/python3.?/site-packages/
   ```
5. Test your installation:
   ```sh
   python multiprot/testing.py -v 2
   ```
   This will run the complete multiprot test suite (can take some minutes). If both Ranch and Pulchra are installed, it should finish without errors.

# Using multiprot on the commandline

Go to the package directory and run the `multiprot_script.py` script with the required arguments.
To see a full list of arguments run `python multiprot_script.py -h`. Here you can find a few examples, where all the file names are relative to the `multiprot/multiprot` package folder:

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
