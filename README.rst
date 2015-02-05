Loop modeling benchmark
=======================
The goal of loop modeling is to predict the native conformation of an internal 
stretch of protein backbone.  This is an important problem in areas such as  
homology modeling, protein design, and structure determination.  This benchmark 
comprises a standard set of well-curated loops for the purpose of comparing 
different loop modeling algorithms.  Each algorithm predicts hundreds of 
structures for each loop and is judged based on the backbone RMSD between those 
predictions and the known structure.

Licensing
---------
The contents of the repository *where possible* are licensed under the Creative 
Commons Attribution 4.0 International License. The license only applies to 
files which either: i) include the license statement; or ii) which are 
explicitly listed in some file in the repository as being covered by the 
license. All other files may be covered under a separate license. The LICENSE 
file in the root of this repository is present only for the convenience of the 
user to indicate the license which covers any novel content presented herein.

Downloading the benchmark
-------------------------
The benchmark is hosted on GitHub. The most recent version can be checked out 
using the `git <http://git-scm.com/>`_ command-line tool:

::

  git clone https://github.com/Kortemme-Lab/loop_modeling.git

Directories in this archive
---------------------------
This archive contains the following directories:

- *libraries* : contains common code used by all the scripts comprising the 
  benchmark.
- *input* : contains the input files for the benchmark.
- *output* : these directories are empty by default. This is the default output 
  location for protocols if they are run on the local machine.
- *output/sample* : contains sample output data that can be used to test the 
  analysis script.
- *analysis* : contains the analysis script used to analyze the output of a 
  prediction run. All protocols are expected to produce output in a format 
  compatible with the analysis script (see analysis/README.rst for details of 
  this format).
- *protocols* : contains the scripts needed to run a prediction for each 
  protocol.
- *hpc* : contains scripts that can be used to run the entire benchmark using 
  specific cluster architectures. For practical reasons, a limited number of 
  cluster systems are supported. Please feel free to provide scripts which run 
  the benchmark for your particular cluster system.

Most of these directories contain README files.
