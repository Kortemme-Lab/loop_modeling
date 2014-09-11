#!/usr/bin/env python2

import sys
import os
import shutil
import subprocess
import glob

from libraries import utilities; utilities.require_chef()
from libraries import settings; settings.load()
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

def run_benchmark(script, pdbs,
        name=None, desc=None, vars=(), fast=False, explicit_log=False):

    pdbs = [x for x in sorted(pdbs)]

    # Make sure all the inputs actually exist.

    for pdb in pdbs:
        if not os.path.exists(pdb):
            raise ValueError("'{0}' does not exist.".format(pdb))

    # Create an entry in the benchmarks table.

    with database.connect() as session:
        benchmark = database.Benchmarks(name, desc)

        for pdb in pdbs:
            benchmark_input = database.BenchmarkInputs(pdb)
            benchmark.input_pdbs.append(benchmark_input)

        session.add(benchmark); session.flush()
        benchmark_id = str(benchmark.id)

    # Make sure the log file directories exist.

    if explicit_log:
        utilities.clear_directory('explicit_log')
        log_args = '-o', 'explicit_log', '-j', 'y'
    else:
        log_args = '-o', '/dev/null', '-j', 'y'
        
    # Submit the benchmark to the cluster.

    vars_args = ()
    for var in vars:
        vars_args += '--var', var

    for pdb in pdbs:
        if fast:
            command = (
                    'qsub', '-t', '1-10', '-l', 'h_rt=0:30:00') + log_args + (
                    'benchmark.py', script, pdb, '--id', benchmark_id, '--fast',) + vars_args
        else:
            command = (
                    'qsub', '-t', '1-500') + log_args + (
                    'benchmark.py', script, pdb, '--id', benchmark_id) + vars_args

        subprocess.call(command)


if __name__ == '__main__':
    import optparse

    # Use optparse because it's available on chef.

    usage = 'launch.py [options] <script> <full|mini|pdbs...|test>'
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('--name', dest='name')
    parser.add_option('--desc', dest='desc')
    parser.add_option('--var', dest='vars', action='append', default=[])
    parser.add_option('--compile-only', '-c', action='store_true', dest='compile_only')
    parser.add_option('--execute-only', '-x', action='store_true', dest='execute_only')
    parser.add_option('--explicit-log', action='store_true', dest='explicit_log')
    options, arguments = parser.parse_args()

    if len(arguments) == 0:
        print 'Usage: ' + usage
        print
        print 'launch.py: error: must specify a RosettaScript to benchmark.'
        sys.exit(1)

    if len(arguments) == 1:
        print 'Usage: ' + usage
        print
        print 'launch.py: error: must specify a set of PDB files to benchmark.'
        sys.exit(1)

    script = arguments[0]
    benchmark = set(arguments[1:])

    # Compile rosetta.

    if not options.execute_only:
        error_code = compile_rosetta()
        if error_code != 0:
            sys.exit(error_code)

    if options.compile_only:
        sys.exit(1)

    # Configure the benchmark run.  If a test run is specified, an easy 
    # structure is used and only a few iterations are run.  If a benchmark run 
    # is specified, any structures specified on the command-line will be used.  
    # The 'full' and 'mini' keywords automatically load common structures.

    if 'test' in arguments:
        #pdbs = ['structures/1srp.pdb', 'structures/1bn8.pdb']
        pdbs = ['structures/1srp.pdb']
        fast = True

    else:
        pdbs = set()
        fast = False

        if 'full' in benchmark:
            benchmark.remove('full')
            pdbs |= set(glob.glob('structures/*.pdb'))

        if 'mini' in benchmark:
            benchmark.remove('mini')
            pdbs |= set([
                    "structures/1c5e.pdb",
                    "structures/1cyo.pdb",
                    "structures/1exm.pdb",
                    "structures/2cpl.pdb",
                    "structures/3cla.pdb",
                    "structures/1oth.pdb"])

        pdbs |= benchmark

    run_benchmark(script, pdbs,
            name=options.name, desc=options.desc,
            vars=options.vars, explicit_log=options.explicit_log, fast=fast)

