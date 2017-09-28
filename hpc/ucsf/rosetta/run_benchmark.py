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
    run_benchmark.py --complete ID [options]
    run_benchmark.py --compile-only

Arguments:
    <name>
        The name for this benchmark.  Several benchmark runs may  
        have the same name. In case of ambiguity, the analysis script will 
        automatically pick the most recent run.

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
        This feature has not been implemented for the disk platform.

    --complete ID -c ID
        If benchmark ID is missing data, this command queues up extra jobs for the benchmark using the same settings.

    --keep-old-data
        When a benchmark run is submitted to the cluster, the stdout and stderr files from previous runs are deleted
        from the output directory. To prevent this deletion from occurring, use this flag.

    --use-database
        Report the data to a database if specified. Otherwise dump the data to the disk directly.

    --use-native-structure
        Pass the native structure to Rosetta.
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
from libraries.dataController import DataController 
from check_progress import get_progress


def compile_rosetta():
    rosetta_path = os.path.abspath(os.path.expanduser(settings.rosetta))

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
            # Use newer version of compilers to support c++11
            # The site.settings file of Rosetta should have /opt/rh/devtoolset-4/root/usr/bin
            # prepended to PATH
            'source /opt/rh/devtoolset-4/enable',
            'nohup nice scl enable python27 ./scons.py bin -j16 mode=release extras=mysql'])

    return subprocess.call(compile_command)


def run_benchmark(name, script, pdbs,
        vars=(), flags=None, fragments=None, nstruct=None,
        desc=None, fast=False, non_random=True, use_database=False,
        complete_run=False, keep_old_data=False, use_native_structure=False):

    if nstruct is not None:
        try: nstruct = int(nstruct)
        except: pass
        assert isinstance(nstruct, int)
    if fast:
        nstruct  = nstruct or 10
    else:
        nstruct  = nstruct or 500

    pdbs = [x.strip() for x in sorted(pdbs) if x.strip()]

    # Make sure all the inputs actually exist.

    for pdb in pdbs:
        if not os.path.exists(pdb):
            raise ValueError("'{0}' does not exist.".format(pdb))

    # Figure out which version of rosetta is being used.

    git_commit = subprocess.check_output(
            shlex.split('git rev-parse HEAD'), cwd=settings.rosetta).strip()
    git_diff = subprocess.check_output(
            shlex.split('git diff'), cwd=settings.rosetta).strip()

    # Create the the benchmark
    data_controller = DataController('database') if use_database else DataController('disk')
    benchmark_define_dict = { 'name':name, 'script':script, 'nstruct':nstruct,
                              'user':getpass.getuser(), 'desc':desc,
                              'vars':vars, 'flags':flags, 'fragments':fragments,
                              'git_commit':git_commit, 'git_diff':git_diff,
                              'fast':fast, 'non_random':non_random }
    benchmark_id = data_controller.create_benchmark(benchmark_define_dict,
                                                    pdbs)
    
    # Submit the benchmark to the cluster.
    submit_benchmark(benchmark_id, nstruct * len(pdbs), fast=fast, use_database=use_database,
                    complete_run=complete_run, keep_old_data=keep_old_data, 
                    use_native_structure=use_native_structure)


def submit_benchmark(benchmark_id, num_tasks, fast=False, use_database=False, complete_run=False,
                    keep_old_data=False, use_native_structure=False):
    '''Submit the benchmark to the cluster.'''

    qsub_command = 'qsub',
    benchmark_command = ('loop_benchmark.py', str(benchmark_id),
                        '--use-database' if use_database else '--not-use-database',
                        '--complete-run' if complete_run else '--not-complete_run',
                        '--use-native-structure' if use_native_structure else '--not-use-native-structure')
    if fast:
        qsub_command += '-t', '1-{0}'.format(num_tasks)
        qsub_command += '-l', 'h_rt=0:30:00'
    else:
        qsub_command += '-t', '1-{0}'.format(num_tasks)
        qsub_command += '-l', 'h_rt=24:00:00'

    if not keep_old_data:
        utilities.clear_directory('job_output')
    qsub_command += '-o', 'job_output', '-e', 'job_output'
    qsub_cwd = os.path.dirname(__file__)

    subprocess.call(qsub_command + benchmark_command, cwd=qsub_cwd)


def resume_benchmark(benchmark_id, nstruct=None, use_database=False):
    
    # This feature has not been implemented on the disk platfrom
    if not use_database:
        raise Exception("resume benchmark has not been implemented on the disk platform!") 

    qsub_command = 'qsub',
    benchmark_command = ('loop_benchmark.py', benchmark_id, 
                        '--use-database' if use_database else '--not-use-database', 
                        '--not-complete-run')

    # You get weird errors if you forget to cast nstruct from string to int.

    if nstruct is not None: nstruct = int(nstruct)

    # Read the job parameters from the saved data

    data_controller = DataController('database') if use_database else DataController('disk')
    benchmark_define_dict = data_controller.get_benchmark_define_dict(benchmark_id)
    num_pdbs = len(benchmark_define_dict['input_pdbs'])

    # Make sure the right version of rosetta is being used.

    git_commit = subprocess.check_output(
            shlex.split('git rev-parse HEAD'),
            cwd=settings.rosetta).strip()

    git_diff = subprocess.check_output(
            shlex.split('git diff'),
            cwd=settings.rosetta).strip()

    if benchmark_define_dict['git_commit'] != git_commit:
        message = "Benchmark \"{0}\" was run with rosetta commit #{1}, but commit #{2} is currently checked out.  Press [Ctrl-C] to abort or [Enter] to continue."
        message = textwrap.fill(message.format(benchmark_define_dict['id'], benchmark_define_dict['git_commit'][:8], git_commit[:8]))
        raw_input(message)

    elif benchmark_define_dict['git_diff'] != git_diff:
        message = "Uncommitted changes have been made to rosetta since benchmark \"{0}\" was run.  Press [Ctrl-C] to abort or [Enter] to continue."
        message = textwrap.fill(message.format(benchmark_define_dict['id']))
        raw_input(message)

    # Build the qsub command.

    if benchmark_define_dict['fast']:
        qsub_command += '-t', '1-{0}'.format((nstruct or 10) * num_pdbs)
        qsub_command += '-l', 'h_rt=24:00:00'
    else:
        qsub_command += '-t', '1-{0}'.format((nstruct or 500) * num_pdbs)
        qsub_command += '-l', 'h_rt=24:00:00'

    print "Your benchmark \"{0}\" (id={1}) is being resumed".format(
            benchmark_define_dict['name'], benchmark_id)

    # Submit the job.

    utilities.clear_directory('job_output')
    qsub_command += '-o', 'job_output', '-e', 'job_output'

    subprocess.call(qsub_command + benchmark_command)


def complete_benchmark(benchmark_id, nstruct=None, use_database=False, use_native_structure=False):
    name = benchmark_id
    # You get weird errors if you forget to cast nstruct from string to int.

    # Get the DataController of the benchmark
    data_controller = DataController('database') if use_database else DataController('disk')
    
    # Get the progress data for the job
    progress_data = get_progress(data_controller, settings.db_name, benchmark_id)

    # Set up nstruct
    nstruct = progress_data['nstruct']
    if not nstruct:
        sys.exit('The nstruct variable is not set for this benchmark. Exiting.')

    if use_database: 
        # Set up the bins for structures that need extra jobs to be run. We run extra jobs in case these fail as well.
        bins = {5 : [], 10 : [], 20 : [], 30 : []}
        d_bins = bins.keys()
        for input_tag, finished_count in progress_data['CountPerStructure'].iteritems():
            if finished_count < nstruct:
                missing_count = nstruct - finished_count
                if missing_count <= 2:
                    bins[5].append(input_tag)
                elif missing_count <= 5:
                    bins[10].append(input_tag)
                elif missing_count <= 10:
                    bins[20].append(input_tag)
                elif missing_count <= 15:
                    bins[30].append(input_tag)
                else:
                    bin_size = ((int((missing_count - 11)/20.0) + 2) * 20) + 10
                    bins[bin_size] = bins.get(bin_size, [])
                    bins[bin_size].append(input_tag)
        for d_bin in d_bins:
            if not bins[d_bin]:
                del bins[d_bin]

        benchmark_variables = data_controller.get_benchmark_variables( benchmark_id )
        for x in range(0, len(benchmark_variables['rosetta_script_vars']) - 1):
            if benchmark_variables['rosetta_script_vars'][x] != benchmark_variables['rosetta_script_vars'][x + 1]:
                sys.exit('Exception (ambiguity): The benchmark {0} has multiple RosettaScript variable values associated with previous runs: "{1}".'.format(benchmark_id, '", "'.join(map(str, sorted(benchmark_variables['rosetta_script_vars'])))))

        for k, v in sorted(benchmark_variables.iteritems()):
            if len(v) == 0:
                sys.exit('Exception (missing data): The benchmark {0} has no {1} values associated with previous runs.'.format(benchmark_id, k.replace('_', ' ')))
            elif k == 'rosetta_script_vars':
                benchmark_variables[k] = benchmark_variables[k][0]
            elif len(v) > 1:
                sys.exit('Exception (ambiguity): The benchmark {0} has multiple {1} values associated with previous runs: "{2}".'.format(benchmark_id, k.replace('_', ' '), '", "'.join(sorted(v))))
            else:
                benchmark_variables[k] = v.pop()

        for nstruct, pdbs in reversed(sorted(bins.iteritems())): # start the longer jobs first
            run_benchmark(name, benchmark_variables['rosetta_script'], pdbs, vars=benchmark_variables['rosetta_script_vars'],
                          flags=benchmark_variables['rosetta_flags'], fragments=benchmark_variables['rosetta_fragments'], nstruct=nstruct,
                          desc=None, fast=benchmark_variables['fast'], non_random=benchmark_variables['non_random'])

    else:
        unfinished_task_list = data_controller.get_unfinished_task_list(progress_data['MostRecentID'],
                                                                        progress_data['TotalCount'])
        data_controller.create_task_completion_list_file(progress_data['MostRecentID'], unfinished_task_list)
        benchmark_define_dict = data_controller.get_benchmark_define_dict(progress_data['MostRecentID'])
        if len(unfinished_task_list) > 0:
            submit_benchmark(progress_data['MostRecentID'], len(unfinished_task_list), fast=benchmark_define_dict['fast'],
                            use_database=use_database, complete_run=True, keep_old_data=True, use_native_structure=use_native_structure)


if __name__ == '__main__':
    try:
        # Parse command-line options.

        from libraries import docopt
        arguments = docopt.docopt(__doc__)
        #utilities.require_chef()
        settings.load()

        # Compile rosetta.

        if arguments['--complete']:
            complete_benchmark(arguments['--complete'], arguments['--use-database'], arguments['--use-native-structure'])
            sys.exit(1)

        if not arguments['--execute-only']:
            error_code = compile_rosetta()
            if error_code != 0:
                sys.exit(error_code)

        if arguments['--compile-only']:
            sys.exit(1)

        # Decide whether to start a new benchmark or to resume an old one.

        if arguments['--resume'] is not None:
            benchmark_id = arguments['--resume']
            resume_benchmark(benchmark_id, arguments['--nstruct'], arguments['--use-database'])
        
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
                    use_database=arguments['--use-database'],
                    keep_old_data=arguments['--keep-old-data'],
                    use_native_structure=arguments['--use-native-structure']
            )

    except KeyboardInterrupt:
        print


