Input data
==========
The structures comprising the loop modeling benchmark were compiled from 
datasets published by Wang et al. and Sellers et al. (references below).  Each 
structure in the loop benchmark is associated with two files: a PDB file and a 
loop file.  The PDB files specify coordinates for entire proteins, including 
the loop regions to be modeled.  These coordinates have been pre-minimized in 
the rosetta force field.

The loop files specify the regions to be modeled.  The protocols are expected 
to discard coordinate information from within these regions and to predict a 
loop conformation based only on the structure of the rest of the protein.  The 
loop files are in the file format rosetta uses to identify loops, although 
unfortunately this format is opaque and poorly documented.  Each loop file 
consists of 5 numbers:  The first is the start of the loop, the second is where 
the loop should be broken (for protocols that need this information), the third 
is the end of the loop, the fourth is the "skip rate" (which only matters when 
multiple loops are being sampled), and the fifth is a boolean indicating 
whether or not the initial loop coordinates should be discarded.  Although this 
format is specific to rosetta, it should be easy to for other protocols to 
parse as well.

The sets of structures constituting benchmark sets are specified in the *.pdbs 
files.  These files are simply newline-delimited lists of PDB files to be used 
in a benchmark.  Their file extension does not matter.  The reason for having 
multiple benchmark sets is mostly for debugging; during development it's nice 
to be able to run a reduced set of structures to get a quick idea for how a 
protocol is performing.

Table of Contents
=================

full.pdbs
    The set of all the structures that make up the canonical loop benchmark.  
    This benchmark was originally published by Mandell et al. in the paper 
    describing the KIC algorithm.  This dataset is a subset of the combined 
    Wang et al. and Sellers et al. datasets.

mini.pdbs
    A reduced set of structures that useful for getting a rough view of the 
    performance of an algorithm.  A mix of easy, medium, and hard structures 
    are included in this set.

ions.pdbs
    The set of all structures containing metal ions.  In most cases the ions 
    are not located near the loop being sampled.  This set maybe useful for 
    debugging a protocol that's having a hard time loading ions.

structures
    The directory containing all the PDB and loop files used by the benchmark.  

preparation
    In principle, the directory containing the scripts that were used to 
    prepare the benchmark structures.  In practice, these scripts don't work 
    right now and I don't know exactly how the structures were prepared.  See 
    ``preparation/README.rst`` for more information.

fragments
    Fragment files used by some of the protocols (e.g. CCD and KIC with 
    fragments)  These files contain coordinates of peptide fragments in the PDB 
    with sequence similarity to the loops in the benchmark.

References
==========
Sellers BD, Zhu K, Zhao S, Friesner RA, Jacobson MP (2008) Toward better 
refinement of comparative models: predicting loops in inexact environments.  
Proteins 72: 959–971. doi: 10.1002/prot.21990

Wang C, Bradley P, Baker D (2007) Protein-protein docking with backbone 
flexibility. Journal of molecular biology 373: 503–519. doi: 
10.1016/j.jmb.2007.07.050 
