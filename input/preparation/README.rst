minpack.sh, preprocess.sh, preprocess_all.sh input preparation scripts
======================================================================

These scripts were used to prepare all of the input files for Rosetta in a uniform way however they are
currently missing some functionality. We have archived them mainly to include the Rosetta flags used in
the preparation. The structures in input/structures/*loop length*/rosetta/reference have been minimized in the Rosetta
all-atom force field (score function).

remove_native_information.py input preparation script
=====================================================

This script was used to remove the loop residues and surrounding chains from both the original
(input/structures/*loop length*/rcsb/reference) and the Rosetta-minimized (input/structures/*loop length*/rosetta/reference) structures. The
resulting structures are stored in the respective sibling directory named "pruned". This was done to remove
possible native conformational bias in the loop region.

To allow the older Rosetta protocols to work with these modified structures, we added the loop residue
backbone atoms (N, CA, and C only) back into the pruned, preminimized Rosetta structures in a non-native
conformation (almost linearly and equidistantly spacing backbone atoms). This code can also
be found in this script and could be adapted for other pruned input.

This script relies on the Kortemme Lab tools repository which will be made at https://github.com/Kortemme-Lab/klab.

