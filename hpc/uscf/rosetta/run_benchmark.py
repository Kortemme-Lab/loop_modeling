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
Launch a loops benchmark run.  Options are provided so you can easily control 
every important aspect of the benchmark, including which protocol to test, 
which structures to use, and how much simulation to do one each one.  The 
benchmark results are written to the MySQL database specified in the settings 
file to facilitate storage and organization.  This script automatically 
compiles rosetta with database support before each run.

Usage:
    run_benchmark.py <name> <script> <pdbs>... [--var=VAR ...] [options]
    run_benchmark.py --resume ID [options]

Arguments:
    <name>
        The name for this benchmark.  It's ok for several benchmark runs to 
        have the same name.  The analysis script will automatically pick the 
        most recent one when given an ambiguous name.

    <script>
        A rosetta XML script to execute.  Commonly used scripts can be found in 
        the benchmarks directory of this repository.  Feel free to add more!

    <pdbs>
        A list of PDB files to use for the benchmark.  Files with the extension 
        '*.pdbs' are expected to contain a list of PDB files to include (one on 
        each line).  Commonly used PDB lists can be found in the benchmarks 
        directory of this repository.

Options:
    --desc DESC -m DESC
        Give a more detailed description of this benchmark run.

    --compile-only
        Compile rosetta but don't run the benchmark.

    --execute-only -x
        Launch the benchmark without compiling rosetta.  I never use this flag 
        when launching full-scale benchmarks, but for test runs it's not worth 
        waiting 2-3 minutes for scons to figure out that nothing has changed.

    --var VAR
        Specify a rosetta-scripts macro substitution to make.  This option can 
        be specified any number of times.  Each instance of this option should 
        specify and name and a value like so: "--var name=value".

    --flags OPT
        Specify a rosetta flag file containing extra options for this run.

    --fragments DIR
        Specify which fragments files to use for this benchmark.  If this flag 
        is not specified, no fragments files will be given to rosetta, which 
        may cause some XML scripts to crash.

    --nstruct NUM -n NUM
        Specify how many simulations to do for each structure in the benchmark.
        The default value is 500.

    --fast
        Run jobs with a very small number of iterations and lower the default 
        value of --nstruct to 10.  This is useful when you're just making sure 
        a new algorithm runs without crashing.

    --non-random
        Use the SGE Task ID as the random seed for each job.  Two jobs both run 
        with this flag should produce the exact same trajectories.

    --resume ID -r ID
        Expand the given benchmark by running more jobs.  The new jobs will use 
        the same PDB files, rosetta script files, rosetta script variables, 
        rosetta flag files, and "fast" settings as the previous jobs did.  
        However, results may differ if the contents of these files, or the 
        checked out version of rosetta, are changed.
"""

import getpass
import glob
import json
import os
import shlex
import shutil
import subprocess
import sys

from libraries import utilities
from libraries import settings
from libraries import database


def compile_rosetta():
    rosetta_path = os.path.abspath(settings.rosetta)

    # Setup the compiler for the cluster.

    qb3_settings = os.path.join(
            rosetta_path, 'source', 'tools', 'build', 'site.settings.qb3')
    site_settings = os.path.join(
            rosetta_path, 'source', 'tools', 'build', 'site.settings')

    shutil.copyfile(qb3_settings, site_settings)

    # Copy the mysql header files into rosetta.

    mysql_headers = (
            '/netapp/home/kbarlow/lib/'
            'mysql-connector-c-6.1.2-linux-glibc2.5-x86_64/include/*')
    rosetta_headers = os.path.join(
            rosetta_path, 'source', 'external', 'dbio', 'mysql')
    
    for source_path in glob.glob(mysql_headers):
        target_name = os.path.basename(source_path)
        target_path = os.path.join(rosetta_headers, target_name)
        if not os.path.exists(target_path):
            os.symlink(source_path, target_path)

    # Compile rosetta.

    scons_path = os.path.join(rosetta_path, 'source')

    compile_command = 'ssh', 'iqint', '; '.join([
            'cd "%s"' % scons_path,
            'nohup nice ./scons.py bin -j16 mode=release extras=mysql'])

    return subprocess.call(compile_command)

def run_benchmark(name, script, pdbs,
        vars=(), flags=None, fragments=None, nstruct=None,
        desc=None, fast=False, non_random=True):

    pdbs = [x for x in sorted(pdbs)]

    # Make sure all the inputs actually exist.

    for pdb in pdbs:
        if not os.path.exists(pdb):
            raise ValueError("'{0}' does not exist.".format(pdb))

    # Figure out which version of rosetta is being used.

    git_commit = subprocess.check_output(
            shlex.split('git rev-parse HEAD'), cwd=settings.rosetta).strip()
    git_diff = subprocess.check_output(
            shlex.split('git diff'), cwd=settings.rosetta).strip()

    # Test the database connection

    try: database.test_connect()
    except RuntimeError, error: 
        print error
        sys.exit(1)
    
    # Create an entry in the benchmarks table.
    with database.connect() as session:
        benchmark = database.Benchmarks(
                name, script,
                user=getpass.getuser(), desc=desc,
                vars=json.dumps(vars), flags=flags, fragments=fragments,
                git_commit=git_commit, git_diff=git_diff,
                fast=fast, non_random=non_random,
        )

        for pdb in pdbs:
            benchmark_input = database.BenchmarkInputs(pdb)
            benchmark.input_pdbs.append(benchmark_input)

        session.add(benchmark); session.flush()
        benchmark_id = str(benchmark.id)

    print "Your benchmark \"{0}\" (id={1}) has been created".format(
            name, benchmark_id)

    # Submit the benchmark to the cluster.

    qsub_command = 'qsub',
    benchmark_command = 'loop_benchmark.py', benchmark_id

    if nstruct is not None:
        assert isinstance(nstruct, int)

    if fast:
        qsub_command += '-t', '1-{0}'.format((nstruct or 10) * len(pdbs))
        qsub_command += '-l', 'h_rt=0:30:00'
    else:
        qsub_command += '-t', '1-{0}'.format((nstruct or 500) * len(pdbs))
        qsub_command += '-l', 'h_rt=6:00:00'

    utilities.clear_directory('job_output')
    qsub_command += '-o', 'job_output', '-e', 'job_output'
    qsub_cwd = os.path.dirname(__file__)

    subprocess.call(qsub_command + benchmark_command, cwd=qsub_cwd)

def resume_benchmark(benchmark_id, nstruct=None):
    qsub_command = 'qsub',
    benchmark_command = 'loop_benchmark.py', benchmark_id

    # You get weird errors if you forget to cast nstruct from string to int.

    if nstruct is not None: nstruct = int(nstruct)

    # Read the job parameters from the database.

    with database.connect() as session:
        benchmark = session.query(database.Benchmarks).get(benchmark_id)
        num_pdbs = len(benchmark.input_pdbs)

        # Make sure the right version of rosetta is being used.

        git_commit = subprocess.check_output(
                shlex.split('git rev-parse HEAD'),
                cwd=settings.rosetta).strip()

        git_diff = subprocess.check_output(
                shlex.split('git diff'),
                cwd=settings.rosetta).strip()

        if benchmark.git_commit != git_commit:
            message = "Benchmark \"{0}\" was run with rosetta commit #{1}, but commit #{2} is currently checked out.  Press [Ctrl-C] to abort or [Enter] to continue."
            message = textwrap.fill(message.format(benchmark.id, benchmark.git_commit[:8], git_commit[:8]))
            raw_input(message)

        elif benchmark.git_diff != git_diff:
            message = "Uncommitted changes have been made to rosetta since benchmark \"{0}\" was run.  Press [Ctrl-C] to abort or [Enter] to continue."
            message = textwrap.fill(message.format(benchmark.id))
            raw_input(message)

        # Build the qsub command.

        if benchmark.fast:
            qsub_command += '-t', '1-{0}'.format((nstruct or 10) * num_pdbs)
            qsub_command += '-l', 'h_rt=0:30:00'
        else:
            qsub_command += '-t', '1-{0}'.format((nstruct or 500) * num_pdbs)
            qsub_command += '-l', 'h_rt=4:00:00'

        print "Your benchmark \"{0}\" (id={1}) is being resumed".format(
                benchmark.name, benchmark_id)

    # Submit the job.

    utilities.clear_directory('job_output')
    qsub_command += '-o', 'job_output', '-e', 'job_output'

    subprocess.call(qsub_command + benchmark_command)

if __name__ == '__main__':
    try:
        # Parse command-line options.

        from libraries import docopt
        arguments = docopt.docopt(__doc__)
        utilities.require_chef()
        settings.load()

        # Compile rosetta.

        if not arguments['--execute-only']:
            error_code = compile_rosetta()
            if error_code != 0:
                sys.exit(error_code)

        if arguments['--compile-only']:
            sys.exit(1)

        # Decide whether to start a new benchmark or to resume an old one.

        if arguments['--resume'] is not None:
            benchmark_id = arguments['--resume']
            resume_benchmark(benchmark_id, arguments['--nstruct'])
        
        else:
            name = arguments['<name>']
            script = arguments['<script>']
            pdb_args, pdbs = set(arguments['<pdbs>']), set()

            # Decide which structures to benchmark.

            for path in pdb_args:
                if path.endswith('.pdb') or path.endswith('.pdb.gz'):
                    pdbs.add(path)

                elif path.endswith('.pdbs'):
                    with open(path) as file:
                        pdbs.update(line.strip() for line in file)

                else:
                    print "Unknown input structure '{0}'.".format(path)
                    sys.exit(1)

            for pdb in pdbs:
                if not os.path.exists(pdb):
                    print "Unknown input structure '{0}'.".format(pdb)
                    sys.exit(1)

            # Run the benchmark.

            run_benchmark(
                    name, script, pdbs,
                    vars=arguments['--var'],
                    flags=arguments['--flags'],
                    fragments=arguments['--fragments'],
                    nstruct=arguments['--nstruct'],
                    desc=arguments['--desc'],
                    fast=arguments['--fast'],
                    non_random=arguments['--non-random'],
            )

    except KeyboardInterrupt:
        print


