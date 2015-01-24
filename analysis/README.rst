====================================
Metrics
====================================

This benchmark uses ...two?... different metrics:

- Median Top5 RMSD; the median loop RMSD across all 45 cases using the top 5 lowest energy structures for each case where each RMSD is taken from the structure which is closest to the native; and
- Median % ≤1Å; the median percentage of sub-Ångström models across all 45 cases.

For more complete explanations of the metrics, the reader is referred to the ...year... publication from ...authors... [1]_

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Median Top5 RMSD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

...more_explanation?_maybe_give_motivation...

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Median % ≤1Å
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

...more_explanation?_maybe_give_motivation...


================
Running analysis
================

~~~~~~~~~~~~~~
Required tools
~~~~~~~~~~~~~~

The analysis script requires ...list_of_dependencies...
The script has been tested using ...e.g.Python2.6.6_and_... .

~~~~~~~~~~~~~~~~~~~~~~~~
Expected format of input
~~~~~~~~~~~~~~~~~~~~~~~~

The analysis script expects data to be in one of two formats: i) a relational database following a particular schema
(see below); or ii) a text file following a particular format (see below).

---------------
Database schema
---------------

...pseudo-definition-of-tables_or_CREATE_TABLE_commands_for_MySQL/some_other_engine_...

----------------
Flat file format
----------------

The flat file should be a tab-separated file. The first line in the file should be a header line:

::

  #PDB	Model	Loop_rmsd	Total_energy	Runtime

where the loop RMSD is measured in Ångströms, the energy is measured in the units of the loop modeling protocol, and the
runtime is measured in seconds. Subsequent lines in the file should be data following this format *e.g.*

::

  1a8d	1	6.4012	-496.689	2310
  1a8d	2	5.82274	-505.773	3444
  1a8d	3	5.01459	-504.018	4635
  ...
  1a8d	500	5.58601	-501.161	5071
  1arb	1	2.51994	-289.462	4684
  ...

~~~~~~~~~~~~~~~~~
Report generation
~~~~~~~~~~~~~~~~~

::

  ...some_command_lines...

This should create a LaTeX report ...explain_what_report_contains...

==========
References
==========

.. [1] ...some_publication...
