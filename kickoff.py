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

def run_benchmark(name, script, pdbs,
        desc=None, vars=(), nstruct=None, debug=False, explicit_log=False):

    pdbs = [x for x in sorted(pdbs)]

    # Make sure all the inputs actually exist.

    for pdb in pdbs:
        if not os.path.exists(pdb):
            raise ValueError("'{0}' does not exist.".format(pdb))

    # Create an entry in the benchmarks table.

    with database.connect() as session:
        import getpass; user = getpass.getuser()
        benchmark = database.Benchmarks(name, user=user, desc=desc)

        for pdb in pdbs:
            benchmark_input = database.BenchmarkInputs(pdb)
            benchmark.input_pdbs.append(benchmark_input)

        session.add(benchmark); session.flush()
        benchmark_id = str(benchmark.id)

    print "Your benchmark \"{0}\" (id={1}) has been created".format(
            name, benchmark_id)

    # Submit the benchmark to the cluster.

    qsub_command = 'qsub',
    benchmark_command = 'loop_benchmark.py', benchmark_id, script

    if debug:
        qsub_command += '-t', '1-{0}'.format((nstruct or 10) * len(pdbs))
        qsub_command += '-l', 'h_rt=0:30:00'
    else:
        qsub_command += '-t', '1-{0}'.format((nstruct or 500) * len(pdbs))
        qsub_command += '-l', 'h_rt=3:00:00'

    if explicit_log:
        utilities.clear_directory('explicit_log')
        qsub_command += '-o', 'explicit_log', '-j', 'y'
    else:
        qsub_command += '-o', '/dev/null', '-j', 'y'
        
    if debug:
        benchmark_command += '--fast',
    for var in vars:
        benchmark_command += '--var', var

    subprocess.call(qsub_command + benchmark_command)


if __name__ == '__main__':
    import optparse

    # Use optparse because it's available on chef.

    usage = 'kickoff.py [options] <name> <script> <pdbs>'
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('--desc', '-m', dest='desc')
    parser.add_option('--var', action='append', default=[], dest='vars')
    parser.add_option('--nstruct', '-n', type=int, dest='nstruct')
    parser.add_option('--compile-only', '-c', action='store_true', dest='compile_only')
    parser.add_option('--execute-only', '-x', action='store_true', dest='execute_only')
    parser.add_option('--explicit-log', action='store_true', dest='explicit_log')
    parser.add_option('--debug', action='store_true', dest='debug')
    options, arguments = parser.parse_args()

    if len(arguments) == 0:
        print 'Usage:', usage
        print
        print 'kickoff.py: error: must provide a name for this benchmark.'

    if len(arguments) == 1:
        print 'Usage:', usage
        print
        print 'kickoff.py: error: must specify a RosettaScript to benchmark.'
        sys.exit(1)

    if len(arguments) == 2:
        print 'Usage:', usage
        print
        print 'kickoff.py: error: must specify a set of PDB files to benchmark.'
        sys.exit(1)

    name = arguments[0]
    script = arguments[1]
    pdb_paths = set(arguments[2:])

    # Compile rosetta.

    if not options.execute_only:
        error_code = compile_rosetta()
        if error_code != 0:
            sys.exit(error_code)

    if options.compile_only:
        sys.exit(1)

    # Decide which structures to benchmark.

    pdbs = set()

    for path in pdb_paths:
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
            desc=options.desc, vars=options.vars, nstruct=options.nstruct,
            debug=options.debug, explicit_log=options.explicit_log)

