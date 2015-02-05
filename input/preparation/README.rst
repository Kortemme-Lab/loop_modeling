Input Preparation Scripts
=========================
These scripts were meant to be used to prepare all of the input files in a 
uniform way.  Unfortunately, I had a hard time getting them to work, so I gave 
up and just used the structures that our lab had previously been using for this 
benchmark.  I couldn't find the scripts used to prepare those structures, so I 
don't really know how uniformly they were prepared.

If you ever need to add more structures to the benchmark and want to prepare 
the structures in a uniform way, maybe start by trying to either debug or 
rewrite these scripts.  If you come up with a good, script-able way to prepare 
the structures, feel free to apply those scripts to all the existing 
structures.

Also note that these preparation scripts are specific to rosetta, especially in 
that the minimize the structures in the rosetta force field.  This is 
appropriate for all the rosetta protocols, but of course is inappropriate for 
any protocol that doesn't use the rosetta score function.  If you are 
interested in adding such a protocol to the benchmark, the structure of the 
input directory may need to change.  The current inputs should be moved into a 
rosetta specific directory, and a new directory should be created for your 
benchmark's pre-minimized inputs.  A new directory containing the raw PDBs, 
downloaded straight from the PDB, should also be created for convenience.

