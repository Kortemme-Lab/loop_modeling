Rosetta Loop Modeling Protocols
===============================
This directory contains rosetta scripts describing a number of different loop 
modeling algorithms.  In general these scripts can do just about anything, but 
for the most part they are variations on the same algorithm:

1. Discard the initial loop coordinates and build new ones that connect the 
   takeoff and landing points.  These coordinate don't need to be very 
   realistic, they just need to be a starting point for refinement.

2. Refine the loop coordinates using the "centroid" score function.  This score 
   function represents all the protein sidechains as spheres and in 
   consequently much smoother than the fullatom score function.

3. Refine the loop coordinates using the "fullatom" score function.  This score 
   function represents every atom, as its name implies, and so is important for 
   predicting atomic-level details.

The variations mostly relate to how the backbone is sampled.  Rosetta is a 
Monte Carlo framework, so different methods use different moves to sample 
backbone conformations.  The most well-represented move in the benchmark is 
KIC, which was described by Mandell et al.  This move is based on an algorithm 
inspired by robotics to analytically calculate closed backbone conformations.  
Some of the protocols also differ in how they ramp the temperature and certain  
score function weights.

Protocol Descriptions
---------------------
kic.xml
    This script carries out "Next Generation KIC" as described by Stein et al.  
    This is based on the analytic KIC move for sampling closed backbones.  Phi 
    and psi torsions are picked from the 2-body Ramachandran distribution, 
    omega torsions are sampled and are picked from a narrow empirical 
    distribution, and both the score function and the temperature are ramped 
    over the course of the simulation.  After each backbone move, the sidechain 
    rotamers are optimized and the entire loop region in minimized by steepest 
    descent before the resulting conformation is subjected to the Metropolis 
    criterion.

kic_with_frags.xml
    This script carries out "KIC with Fragments", which has not yet been 
    published.  It uses a library of fragments culled from the PDB to sample 
    phi, psi, and omega torsions.  Note that you cannot use this script unless 
    you first generate a fragment library.  In the future we plan to include 
    the fragment library in the benchmark.

kic_no_sfxn_ramp.xml
    This is dummy protocol we use just for comparison to other things.

legacy/loopmodel.xml
    This script is a thin wrapper around the legacy rosetta ``loopmodel`` app, 
    which cannot be called directly by this benchmark because it is not a 
    rosetta script.  The behavior of this script is controlled by the *.flag 
    files included in this directory, which can be passed to 
    ``hpc/ucsf/rosetta/run_benchmark.py``.
    
    - ccd.flags: The algorithm described by Wang et al.
    - kic.flags: The original KIC algorithm described by Mandell et al.
    - next_gen_kic.flags: Identical but less flexible than ../kic.xml
    - refactored_kic.flags: Identical but less flexible than ../kic.xml
    - kic_with_frags.flags: Identical but less flexible than ../kic_with_frags.xml

All of the protocols listed above are described more fully in the `loop 
modeling section 
<https://www.rosettacommons.org/docs/latest/loop-modeling-movers.html>`_ of the 
rosetta documentation.

Debugging Protocols
-------------------
The ``debug_protocol.py`` script is provided to make debugging new protocols 
easier.  This script can run individual protocols locally.  By default it runs 
a very limited number of cycles.  It also provides an option to automatically 
drop into the gdb, which is very convenient.  This script should not be used to 
actually collect benchmark data.

References
----------
Mandell DJ, Coutsias EA, Kortemme T. Sub-angstrom accuracy in protein loop 
reconstruction by robotics-inspired conformational sampling. Nature methods 
2009;6(8):551-552. doi:10.1038/nmeth0809-551.

Stein A, Kortemme T. Improvements to robotics-inspired conformational sampling 
in rosetta. PLoS One. 2013 May 21;8(5):e63090. doi: 
10.1371/journal.pone.0063090.

Wang C, Bradley P, Baker D (2007). Protein-protein docking with backbone 
flexibility. J. Mol. Biol. 373, 503-19.
