Input data
==========
The structures comprising the loop modeling benchmark were compiled from 
datasets published by Fiser et al. and Zhu et al. (references below). Each
case in the loop benchmark is associated with two types of input: a PDB file
describing the protein structure and a loop definition identifying the region
where the loop should be inserted.

The structural files (input/structures/rcsb/reference) identified by the datasets
contain the coordinates for the loop region which allows us to determine how closely
a computational method can predict the loop structure. Depending on the method,
including the loop and surrounding sidechain coordinates may bias the prediction
by including native structural information. For this reason, we have included
input files with the loop residues and surrounding sidechains removed (input/structures/rcsb/pruned).
For the sidechain removal, we chose to remove any sidechain with a heavy atom within
10 angstroms of any heavy atom of a loop residue. The code used to remove these
coordinates can be found in the input/preparation directory. This allows us to benchmark
methods using the type of input to which loop modeling methods will commonly be applied
in practice.

Loop definitions specify the regions to be modeled. We have included those in JSON_
format (input/structures/loop_definitions.json). There is one loop definition per
PDB structure and the definition identifies the chain, first residue, and last residue
of the loop where the residue identifiers use the PDB format (columns 23-27 of the coordinate
records). The loop definitions are taken from the Mandell et al. paper and include loops
not including in the benchmark proper as they were removed by Mandell et al. for not
passing criteria defined in that paper (the PassedFilter field indicates whether the loop
should be considered in this benchmark). The definition also includes the canonical amino acid
sequence and a reference to the publication which included the case in its dataset.

We also provide input files for use with Rosetta methods which have been preminimized in the
Rosetta all-atom force field (input/structures/rosetta). First, the full structures were
minimized in the force field / score function (input/structures/rosetta/reference). The
Rosetta score terms are included for each structure. Next, we removed the loop residues from
the structures as we did for the original structures (input/structures/rosetta/pruned). The
Rosetta methods considered here perform optimally if the N and CA atoms of the first loop
residue and the CA and C atoms of the last loop residue are retained so we have done so. The
information needed to rebuild the loop (the loop definition and sequence) are provided in a
JSON and FASTA file for each case, with the JSON files again using PDB residue identifiers.
Finally, a set of structures for use with particular Rosetta loop modeling methods (CCD, KIC, NGK)
are provided (input/structures/rosetta/kic). We added the loop residue backbone atoms back
into these structures but in an arbitrary and non-native conformation and provided Rosetta
loop definition files using Rosetta residue numbering (consecutive positive, non-zero integers).
These structures are currently used to benchmark certain Rosetta methods.

The Rosetta loops_ file format is particular to that software suite. It specifies the start and
end loop residues using Rosetta residue numbering and can additionally specify a cut point,
a "skip rate" (this only matters when multiple loops are being sampled), and a boolean indicating
whether or not the initial loop coordinates should be discarded.

The sets of structures constituting benchmark sets are specified in the \*.pdbs files. These
files are simply newline-delimited lists of PDB files to be used in a benchmark. Their
file extension does not matter. The reason for having multiple benchmark sets is mainly for
the purposes of debugging and development.

Notes
=====

Seven of the benchmark PDB files are either homodimers or homotrimers. We have kept all chains
in the included files. When coordinates were removed, they were only removed from one chain (always
chain A). Some non-canonical residues in the RCSB structures (*e.g.* PCA) were removed in the Rosetta
structures as were any residues with empty occupancy (*e.g.* in 2pia). There were a limited number of
atoms with coordinates for multiple conformations. For those atoms, the Rosetta structures use coordinates
from the A conformation (the occupancies are typically close to 0.5 in these cases).


Table of Contents
=================

full.pdbs
    The set of all the structures that make up the canonical loop benchmark.  
    This benchmark was originally published by Mandell et al. in the paper 
    describing the KIC algorithm.  This dataset is a subset of the combined 
    Fiser et al. and Zhu et al. datasets.

mini.pdbs
    A reduced set of structures that useful for getting a rough view of the 
    performance of an algorithm.  A mix of easy, medium, and hard structures 
    are included in this set.

ions.pdbs
    The set of all structures containing metal ions.  In most cases the ions 
    are not located near the loop being sampled.  This set maybe useful for 
    debugging a protocol that's having a hard time loading ions.

full\_*.pdbs, mini\_*.pdbs ions\_*.pdbs
    Versions of the above lists specific to particular methods.

structures
    The directory containing all the PDB and loop files used by the benchmark.  
    The .pdb files are in the standard (RCSB) numbering.

structures/loop_definitions.json
    Contains the definition of all loops as given in Mandell, Coutsias, &
    Kortemme (doi:10.1038/nmeth0809-551) in PDB numbering.

structures/rcsb/reference
    The PDB files as downloaded from the RCSB website.

structures/rcsb/pruned
    The RCSB files with the loop residues and surrounding sidechains removed. These
    are the input files for generic methods.

structures/rosetta/reference
    The original RCSB files minimized in the Rosetta force field.

structures/rosetta/pruned
    The preminimized structures above with the loop residues and surrounding sidechains removed.
    The .loop.json files contain the loop definitions in PDB numbering using a JSON format recognized by Rosetta.

structures/rosetta/kic
    The pruned structures above with the loop residue backbone atoms (N, CA, C only) added in a non-native
    conformation (see ``preparation/README.rst``). These structures are used for the CCD, KIC, and NGK Rosetta
    loop modeling methods. The .loop files contain the loop definitions in the older Rosetta loop file
    format which uses Rosetta numbering.

preparation
    This directory containing the scripts that were used to prepare the benchmark structures. At present, details for
    the preminimization Rosetta step are incomplete. The "pruning" scripts are available here.

fragments
    Fragment files used by some of the protocols (e.g. CCD and KIC with fragments)  These files contain coordinates of
    peptide fragments in the PDB with sequence similarity to the loops in the benchmark.


Adding structures
=================

We provide input structures which have been minimized using a Rosetta score function. Depending
on your method, you may be performing a similar preminimization step before loop building. The
scripts in the preparation directory may be useful in this regard.

If you would like to contribute updates to the preparation scripts or the set of input structures,
please feel free to contact us.


References
==========

Fiser A, Do RK, and Sali A (2000). Modeling of loops in protein structures.
Protein Science 2000 9(9): 1753–1773. doi: 10.1110/ps.9.9.1753
http://salilab.org/decoys/

Mandell DJ, Coutsias EA, Kortemme T (2009). Sub-angstrom accuracy in protein loop
reconstruction by robotics-inspired conformational sampling. Nature methods
2009;6(8):551-552. doi:10.1038/nmeth0809-551.

Sellers BD, Zhu K, Zhao S, Friesner RA, Jacobson MP (2008). Toward better
refinement of comparative models: predicting loops in inexact environments.  
Proteins 72: 959–971. doi: 10.1002/prot.21990
http://www.jacobsonlab.org/decoy.htm

Wang C, Bradley P, Baker D (2007). Protein-protein docking with backbone
flexibility. Journal of molecular biology 373: 503–519. doi: 
10.1016/j.jmb.2007.07.050 

Zhu K, Pincus, DL, Zhao S, Friesner RA (2006). Long loop prediction using the
protein local optimization program. Proteins 65: 438–452. doi: 10.1002/prot.21040

.. _JSON: http://www.json.org

.. _loops: https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/loops-file
