# multiprot
Model proteins by connecting two or more structured domains with disordered linkers. Supports single or multiple-chain models and symmetic structures.

This project wraps and adds functionality to two different programs:
- Ranch from the [Ensemble Optimization Method](https://www.embl-hamburg.de/biosaxs/eom.html)
- [Pulchra](http://www.pirx.com/pulchra/index.shtml)

## multiprot installation

The Multiprot installation requires three steps: (i) and (ii) the installation of the two helper applications (Ranch and Pulchra) and then (iii) the setup of the python package. Here are the detailed instructions:

### Install Ranch helper application

Ranch is part of the ATSAS software package from the EMBL Hamburg (BioSax) team. It was originally a stand-alone program which has now been merged into the [EOM package](https://www.embl-hamburg.de/biosaxs/eom.html). 

1. Create a user account at [ATSAS website account registration](https://www.embl-hamburg.de/biosaxs/atsas-online/register.php)
2. Download ATSAS package for your architecture from [ATSAS download area](https://www.embl-hamburg.de/biosaxs/atsas-online/download.php)
3. Install package (instructions [here](https://www.embl-hamburg.de/biosaxs/manuals/install.html)). E.g. for Debian or Ubuntu:
    ```sh
    sudo dpkg --install ./ATSAS-2.8.4-1_amd64.deb
    ```
    This will, among other tools, put the `ranch` command into your search path. Open a terminal and verify that the program is correctly installed. Typing the `ranch` command should open the interactive ranch prompt:
    ```sh
    ~> ranch
    *******  ------------------------------------------------------  *******
    *******     RANCH - version 2.2 - (r10552)                    ********
    ...
    ```

### Install Pulchra helper application

[Pulchra](http://www.pirx.com/pulchra/index.shtml) is a tool developed by Piotr Rotkiewicz. We use it to add backbone and side chain atoms to the models produced by Ranch (which is only generating C-alpha traces).

1. Download source code into a temporary directory. On Linux sytems, open a new terminal window and type the following commands:
   ```sh
   cd /tmp
   wget  http://www.pirx.com/downloads/pulchra_306.tgz
   tar xvfz pulchra_306.tgz
   cd pulchra_306
   cc -O3 -o pulchra pulchra.c pulchra_data.c -lm
   ```
   (See pulchra README file for detailed instructions). This generates an executable `pulchra` in your current directory. Test it:
   ```sh
   ./pulchra
   ```
   ...should give you the pulchra help screen with all the options for running the program.
2. Move binary file to system-wide search path, so you can call it from any location. On Linux sytems, now type the following:
   ```sh
   sudo mv pulchra /usr/local/bin/
   ```
3. Then test that pulchra can be called system-wide with the next commands:
   ```sh
   cd ~
   pulchra
   ```
   ...should give you the same help screen.
4. Clean up. You can now remove the source code alltogether or move it to a more appropriate location. For example:
   ```sh
   cd /tmp
   rm -r pulchra_306
   rm pulchra_306.tgz
   ```
   
Note, you can also leave the `pulchra` executable in any other folder instead of the location of step 2, and configure multiprot to find it there. 


### Install multiprot Python package

##### Installing Python 3
You need to have an installation of Python 3 to use multiprot. You can download it from [Python's official webpage](https://www.python.org/downloads/) or through a [package manager](https://docs.python-guide.org/starting/installation/). If you're unsure wether python 3 is installed on your system, type `python3 --version`.

##### Installing multiprot

1. multiprot depends on the python 3 version of [biskit](https://github.com/graik/biskit) which is not yet available on pypi so we have to install it by hand. In your command line, go to a directory where you would like to save multiprot's and biskit's source code and package contents, and install as follows:
   - OPTIONAL: Set up and start a virtual environment. Skip if you prefer installing multiprot system-wide.
     ```sh
     virtualenv --python=python3 venv
     source venv/bin/activate
     pip install --upgrade pip
     ```
   Download the latest version of multiprot and biskit, and install:
   ```sh
   git clone https://github.com/graik/biskit.git biskit
   pip install -r biskit/requirements.txt
   pip install ./biskit
   ```

   ```sh
   git clone https://github.com/strubelab/multiprot.git
   pip install ./multiprot
   ```
2. Test your installation:
   ```sh
   python3 multiprot/multiprot/testing.py -v 2
   ```
   This will run the complete multiprot test suite (can take some minutes). If both Ranch and Pulchra are installed, it should finish without errors.

## Using multiprot on the commandline

You can now call from any location on your command line the `multipr` script with the required arguments.
To see a full list of arguments and related help simply run `multipr -h`. Below you can find a few examples, where all the file names for the arguments are relative to the `multiprot/multiprot` package folder:

### Single chain examples
##### [Example 1](examples/example1/example1.png)
Build a simple single-chain protein joining two structured domains (2z6o.pdb and histone.pdb) with a disordered linker (TGTGTGTGTGTGTGTGTGTG), and saving the resulting model to the [/examples/example1/](/examples/example1/) directory
```
$ multipr --chain testdata/2z6o.pdb TGTGTGTGTGTGTGTGTGTG testdata/histone.pdb --destination ../examples/example1
```

##### [Example 2](examples/example2/example2.png)
Build model joining domAB1.pdb (chain A) with domAB1.pdb (chain B) through a Gly linker and saving the resulting models to [examples/example2/](examples/example2)
```
$ multipr --chain testdata/domAB1.pdb:A TGTGTGTGTGTGTGTGTGTG testdata/domAB2.pdb:B --destination ../examples/example2
```

##### [Example 3](examples/example3/example3.png)
Symmetric protein taking a dimer as a symmetry core, and adding a linker and another dimeric domain (four chains total). This might take a bit longer to run.
```
multipr --chain testdata/domAB1.pdb TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB2.pdb:A --symmetry p2 --symtemplate testdata/domAB1.pdb --destination ../examples/example3
```

### Multiple chain examples

##### [Example 4](examples/example4/example4.png)
Two chains bound at three different places by three dimeric interfaces. You have to provide the PDB file for each dimeric domain, with their coordinates previously fixed. For this, you need to open the pdbs in a molecular visualization software such as pymol, and move their coordinates to the positions where you think they should be located with respect to the other domains, so the models can be constructed. You can open the pdb files used for this example in [multiprot/testdata](multiprot/testdata), see the location of the domains with respect to each other, and how they compare to the [final model](examples/example5/).
```
multipr --chain testdata/domAB1.pdb:A TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/2qud_mod.pdb:A TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB2.pdb:A --chain testdata/domAB1.pdb:B TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/2qud_mod.pdb:B TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB2.pdb:B --fixed testdata/domAB1.pdb testdata/domAB2.pdb testdata/2qud_mod.pdb --destination ../examples/example5
```
You will see a warning from Ranch like this:
```
--------------------- ATTENTION! ---------------------
 - The domains specified as fixed may be too far away -
 - to be connected by the linker extracted from the   -
 - sequence. EOM 2.1 may not be able to make the pool.-
 ------------------------------------------------------
 ```
 If the linkers are long enough to bind the fixed domains, the program should be able to make the models. However, it might take a long time (5-10 minutes).

##### [Example 5](examples/example5/example5.png)
Three chains bound to one another. Since this time there is only one binding point between any given two chains, there is no need to fix any of the domains.
```
multipr --chain testdata/2z6o.pdb TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB1.pdb:A TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/histone.pdb --chain testdata/domAB1.pdb:B TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB2.pdb:A --chain testdata/1it2_A.pdb TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB2.pdb:B TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/histone.pdb --destination ../examples/example6
```

##### [Example 6](examples/example6/example6.png)
Four chains different to each other (no symmetry), joined by a tetrameric core.
```
multipr --chain testdata/5agc.pdb:A TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/2z6o.pdb --chain testdata/5agc.pdb:B TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/histone.pdb --chain testdata/5agc.pdb:C TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/1it2_A.pdb --chain testdata/5agc.pdb:D TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/1it2_A.pdb --destination ../examples/example7
```

##### [Example 7](examples/example7/example7.png)
Two chains with three-fold symmetry.
```
multipr --chain testdata/2ei4_mod.pdb TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB1.pdb:A --chain testdata/histone.pdb TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/domAB1.pdb:B TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG testdata/2z6o.pdb --symmetry p3 --symtemplate testdata/2ei4_mod.pdb --destination ../examples/example8
```
