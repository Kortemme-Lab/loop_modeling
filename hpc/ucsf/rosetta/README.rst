Running the Loop Modeling Benchmark
===================================

Description
-----------
The purpose of the loop modeling benchmark is to compare how well different 
protocols can predict the conformation of internal protein loops.  The primary 
benchmark comprises 45 different curated structures, each with one 12-residue 
loop to predict.  Each protocol makes 500 predictions for each loop and is 
judged based on the distribution of those predictions.  For more information on 
the metrics used to compare protocols, consult the ``analysis`` directory.

The rosetta-based protocols are expressed as XML scripts in the ``protocols`` 
directory.  In general these scripts can do just about anything, but for the 
most part they are variations on the same algorithm:

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

More information on all the protocols that are included in this benchmark can 
found in the ``protocols`` directory.

Table of Contents
-----------------
There are three scripts in this directory:

run_benchmark.py
  This is the script that launches the benchmark.  It takes as input a name, a 
  protocol (as a rosetta script), and a list of PDB files.  The name is just 
  used to refer back to the benchmark later and can be pretty much anything.  
  The rosetta script format is described `here 
  <https://www.rosettacommons.org/docs/latest/RosettaScripts.html>`_.  The only 
  requirement for the scripts is that they contain a ReportToDB mover 
  configured exactly as shown here::

    <ReportToDB name="reporter" task_operations="loop">
      <feature name="TotalScoreFeatures" scorefxn="talaris2013"/>
      <feature name="ProteinRMSDNoSuperpositionFeatures"/>
      <feature name="RuntimeFeatures"/>
    </ReportToDB>
    
  Unfortunately this requirement isn't checked automatically, so you have to 
  make sure you've done this right if you're adding a new protocol.  All of the 
  protocols distributed in the benchmark are properly configured.  If you make 
  a mistake here, the benchmark will run but the analysis script will choke.

  The first time you run this script, you will be prompted for a number of 
  settings required by the benchmark.  These include the path to rosetta, the 
  URL to the database used for IO (see next paragraph), and other miscellaneous 
  things.  These settings are then saved in a file called ``settings.conf`` in 
  the root directory of the repository.  This is a regular text file and can be 
  edited by hand if you ever need to change any of your settings.  These is 
  also support for adding new sections to the settings file that override the 
  default settings you were prompted for.  For more information about this, 
  consult ``libraries/settings.py``.
  
  The benchmark is meant to use a MySQL database for IO.  This script fills in 
  a table with the parameters for each job that needs to run, compiles rosetta 
  with support for MySQL (not a trivial thing to do), does a bunch of error 
  checking, and finally submits the benchmark job.  It also has the ability to 
  resume a previous benchmark, which can be really useful when something goes 
  wrong.  Note that this script is written specifically for the cluster at UCSF 
  and will not work anywhere else.  If you are trying to adapt the loop 
  modeling benchmark to your system, you can either try to adapt this script to 
  your needs, to use this script as inspiration, or to roll something new all 
  on your own.

update_rosetta.py
  This script is just a convenient way to update your checkout of rosetta to 
  the most recent version.

loop_benchmark.py
    This is the script that is actually submitted to the cluster to run the 
    benchmark jobs.  It reads its parameters from the MySQL database and runs 
    rosetta.  If you are trying to adapt the benchmark to a different system, 
    this is where you should look to figure out the precise rosetta command 
    lines that are being called.  Note that this script should only be invoked 
    by the cluster scheduler; you are never meant to invoke it yourself.

Example Command Lines
---------------------
Below are some example command-lines showing how you'd launch a few different 
benchmarks::

  # Run the full KIC benchmark:
  ./run_benchmark.py B1 benchmarks/kic.xml input/full.pdbs

  # Run a modified version of the KIC protocol (with score function ramping 
  # disabled) on the full benchmark:
  ./run_benchmark.py B2 benchmarks/kic_no_sfxn_ramp.xml input/full.pdbs

  # Run a limited number of jobs on a small subset of the full benchmark.  This 
  # is useful for debugging protocols:
  ./run_benchmark.py B3 benchmarks/kic.xml input/mini.pdbs --fast

  # Run 50 extra simulations for each structure in the previous run:
  ./run_benchmark.py --resume B3 --nstruct 50

  # Generate a report comparing all these benchmarks:
  cd ../../../analysis
  ./make_report.py B1 B2 B3

References
----------
Mandell DJ, Coutsias EA, Kortemme T. Sub-angstrom accuracy in protein loop 
reconstruction by robotics-inspired conformational sampling. Nature methods 
2009;6(8):551-552. doi:10.1038/nmeth0809-551.
