Loop modeling benchmark
=======================
The goal of loop modeling is to predict the native conformation of an internal 
stretch of protein backbone.  This is an important problem in areas such as  
homology modeling, protein design, and structure determination.  This benchmark 
comprises a standard set of well-curated loops for the purpose of comparing 
different loop modeling algorithms.  Each algorithm predicts hundreds of 
structures for each loop and is judged based on the backbone RMSD between those 
predictions and the known structure.

Downloading the benchmark
-------------------------
The benchmark is hosted on GitHub. The most recent version can be checked out 
using the `git <http://git-scm.com/>`_ command-line tool::

  git clone https://github.com/Kortemme-Lab/loop_modeling.git

Running the benchmark
---------------------
The benchmark is only designed to run without modification on the QB3 cluster 
at UCSF.  On that cluster, the commands to run the benchmark will look 
something like the examples below.  More information on what these commands do 
and how they can be configured is given in the README.rst files in their 
respective directories::

  cd hpc/ucsf/rosetta
  ./run_benchmark B1 benchmarks/kic.xml benchmarks/full.pdbs

  cd ../../../analysis
  ./make_report B1

Directories in this archive
---------------------------
This archive contains the following directories:

libraries
  Contains common code used by all the scripts comprising the benchmark.

input
    Contains the input files for the benchmark.

output
    These directories are empty by default. This is the default output location 
    for protocols if they are run on the local machine.

output/sample
    Contains sample output data that can be used to test the analysis script.

analysis
    Contains the analysis script used to analyze the output of a prediction 
    run.  All protocols are expected to produce output in a format compatible 
    with the analysis script.

protocols
    Contains the scripts needed to run a prediction for each protocol.

hpc
    Contains scripts that can be used to run the entire benchmark using 
    specific cluster architectures. For practical reasons, a limited number of 
    cluster systems are supported. Please feel free to provide scripts which 
    run the benchmark for your particular cluster system.

Licensing
---------
The contents of the repository *where possible* are licensed under the Creative 
Commons Attribution 4.0 International License. The license only applies to 
files which either: i) include the license statement; or ii) which are 
explicitly listed in some file in the repository as being covered by the 
license. All other files may be covered under a separate license. The LICENSE 
file in the root of this repository is present only for the convenience of the 
user to indicate the license which covers any novel content presented herein.

