See Sketch.png for a diagram of the structure to model

See chainX/runX.txt, chainXY/runXY.txt, chainXYZ/runXYZ.txt and *_annotated.fasta files for details

Sequence of actions to model entire molecule:

Model chain X
    1) ranch seqX.fasta -q=10 -x=histone.pdb -x=dom1_AB_mod.pdb -x=histone.pdb -f=no -f=no -f=no -o=no -o=no -o=no --filesuffix=chainX  -w=chainX -i
    2) Move dom1_B out of dom1_A as a separate chain -> 00001chainX_mod.pdb
    3) python fix_residues_v4c.py -> 00001chainX_mod_fixed.pdb

Model chains XY
    1) Move sequence of chain X into dom1_B -> 00001chainX_mod_fixed_mod.pdb
    2) ranch seqXY.fasta -q=10 -x=00001chainX_mod_fixed_mod.pdb -x=dom1_AB_mod.pdb -f=no -f=no -o=no -o=no --filesuffix=chainXY -w=chainXY -i
    3) Move chain X out of dom1_B and dom2_B out of dom2_A -> 00001chainXY_mod.pdb
    4) python fix_residues_v4c.py -> 00001chainXY_mod_fixed.pdb

Model chains XYZ
    1) Move the sequence of chains X and Y into dom2_B -> 00002chainXY_mod_fixed_mod.pdb
    2) ranch seqXYZ.fasta -q=10 -x=histone.pdb -x=00002chainXY_mod_fixed_mod.pdb -x=histone.pdb -f=no -f=no -f=no -o=no -o=no -o=no --filesuffix=chainXYZ -w=chainXYZ -i
    3) Move chain X and Y out of dom2_B -> 00001chainXYZ_mod.pdb
    4) python fix_residues_v4c.py -> 00001chainXYZ_mod_fixed.pdb

FINAL MODEL -> 00001chainXYZ_mod_fixed.pdb