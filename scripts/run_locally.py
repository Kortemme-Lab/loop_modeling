#!/usr/bin/env python2

"""\
Run a loop modeling rosetta script locally.  This is useful for debugging.  
Presently the script assumes that CMake/Ninja was used to compile rosetta, 
because of the name of the executable it looks for, but this assumption could 
be relaxed.  Output is written to an SQLite database instead of a MySQL one.

Usage:
    run_locally.py <script> [<pdb>] [--var=VAR ...] [options]

Options:
    --var=VAR       Rosetta scripts macro substitution: "--var name=value"
    --verbose -v    Print out the rosetta command-line.
    --unabridged    Run the full number of cycles.
"""

import os, docopt, subprocess
from libraries import settings

arguments = docopt.docopt(__doc__)
script_path = arguments['<script>']
pdb_path = arguments['<pdb>'] or 'structures/1srp.pdb'
loop_path = os.path.splitext(pdb_path)[0] + '.loop'

settings.load()

rosetta_path = os.path.abspath(settings.rosetta)
rosetta_scripts = os.path.join(rosetta_path, 'source', 'bin', 'rosetta_scripts')
rosetta_database = os.path.join(rosetta_path, 'database')

rosetta_command = [
        rosetta_scripts,
        '-database', rosetta_database,
        '-in:file:s', pdb_path,
        '-in:file:native', pdb_path,
        '-inout:dbms:database_name', 'loopmodel.db3',
        '-parser:protocol', script_path,
        '-parser:script_vars',
            'loop_file={0}'.format(loop_path),
            'fast={0}'.format('no' if arguments['--unabridged'] else 'yes'),
]         + arguments['--var']

if arguments['--verbose']:
    print '$ ' + ' '.join(rosetta_command)

subprocess.call(rosetta_command)

