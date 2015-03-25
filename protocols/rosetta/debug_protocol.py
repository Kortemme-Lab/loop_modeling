#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Kale Kundert
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""\
Run a loop modeling rosetta script locally.  This is useful for debugging.  
Presently the script assumes that CMake/Ninja was used to compile rosetta, 
because of the name of the executable it looks for, but this assumption could 
be relaxed.  Output is written to an SQLite database instead of a MySQL one.

Usage:
    run_locally.py <script> [<pdb>] [--var=VAR ...] [options]

Arguments:
    <script>
        A rosetta XML script to execute.  Commonly used scripts can be found 
        in the benchmarks directory of this repository.

    <pdb>
        A single PDB file to model.  The loop file is assumed to have the same 
        base name with a '.loop' extension.

Options:
    --var VAR
        Specify a rosetta-scripts macro substitution to make.  This option can 
        be specified any number of times.  Each instance of this option should 
        specify and name and a value like so: "--var name=value".

    --flags OPT
        Specify a rosetta flag file containing extra options for this run.

    --fragments DIR
        Specify where to look for fragments files.  Note that some XML scripts
        require this option and will crash without it.

    --unabridged -u
        Run the full number of cycles.  By default, only a small number of 
        "test cycles" are run.

    --non-random SEED
        Use the given random seed for this job.

    --output DIR -o DIR
        Specify the directory where the job should be run.  The default is the 
        current working directory.

    --debug
        Run the job in gdb.  This assumes that the jobs has been previously 
        compiled in debug mode.
        
    --verbose -v
        Print out the rosetta command-line.

{settings.rosetta_args}

"""

import os, re, glob, subprocess
from libraries import settings, docopt

arguments = docopt.docopt(__doc__.format(**locals()))
script_path = os.path.abspath(arguments['<script>'])
pdb_path = os.path.abspath(arguments['<pdb>'] or 'input/structures/1srp.pdb')
pdb_tag = os.path.splitext(os.path.basename(pdb_path))[0]
loop_path = re.sub('\.pdb(\.gz)?$', '.loop', pdb_path)
flags_path = arguments['--flags']
fragments_path = arguments['--fragments']
output_dir = arguments['--output'] or 'sandbox'

settings.load(arguments)

rosetta_path = settings.rosetta
if os.path.isabs(rosetta_path): pass
elif os.path.isabs(os.path.expanduser(rosetta_path)): # e.g. ~/some/path
    rosetta_path = os.path.expanduser(rosetta_path)
else:
    rosetta_path = os.path.abspath(rosetta_path)
rosetta_path = os.path.normpath(rosetta_path)

rosetta_scripts = os.path.join(rosetta_path, 'source', 'bin', 'rosetta_scripts')
rosetta_database = os.path.join(rosetta_path, 'database')

if not os.path.exists(rosetta_scripts):
    raise Exception('The RosettaScripts binary %(rosetta_scripts)s does not exist. By default, Rosetta executables are named after the operating system and build type e.g. rosetta_scripts.linuxgccrelease. Either: i) rename an appropriate binary to %(rosetta_scripts)s; or ii) use a symbolic link.' % locals())

rosetta_command = [
        rosetta_scripts,
        '-database', rosetta_database,
        '-in:file:s', pdb_path,
        '-in:file:native', pdb_path,
        '-inout:dbms:database_name', 'loopmodel.db3',
        '-overwrite',
        '-parser:protocol', script_path,
        '-parser:script_vars',
            'loop_file={0}'.format(loop_path),
            'fast={0}'.format('no' if arguments['--unabridged'] else 'yes'),
]         + arguments['--var']

if flags_path is not None:
    rosetta_command += ['@', os.path.abspath(flags_path)]

if fragments_path is not None:
    frag_file = lambda tag, size: glob.glob(os.path.abspath(os.path.join(
            fragments_path, '{0}?'.format(tag),
            '{0}?.200.{1}backup.rewrite.gz'.format(tag, size))))[0]

    rosetta_command += [
            '-loops:frag_sizes', '9', '3', '1',
            '-loops:frag_files',
                frag_file(pdb_tag, 9),
                frag_file(pdb_tag, 3),
                'none'
    ]

if arguments['--non-random'] is not None:
    rosetta_command += ['-run:constant_seed', '-run:jran', arguments['--non-random']]

if arguments['--debug']:
    rosetta_command = ['gdb', '--args'] + rosetta_command

if arguments['--verbose']:
    print '$',
    for argument in rosetta_command:
        if argument[0] in '-@': print '\n   ',
        print argument,
    print; print

if not os.path.exists(output_dir): os.mkdir(output_dir)
subprocess.call(rosetta_command, cwd=output_dir)

