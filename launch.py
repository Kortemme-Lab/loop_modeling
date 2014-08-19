#!/usr/bin/env python2

import sys
import os
import shutil
import subprocess
import glob
import conf

def check_hostname():
    from socket import gethostname

    if gethostname() != 'chef.compbio.ucsf.edu':
        raise SystemExit("This script must be run on chef.")

def compile_rosetta():
    rosetta_path = os.path.abspath(conf.rosetta)

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
            'nohup nice ./scons.py bin -j16 mode=debug extras=mysql'])

    return subprocess.call(compile_command)

def run_benchmark(pdbs, debug=False):
    pdbs = [x for x in sorted(pdbs)]

    # Make sure all the inputs actually exist.

    for pdb in pdbs:
        if not os.path.exists(pdb):
            raise ValueError("'{0}' does not exist.".format(pdb))

    # Wipe the output directory.

    if os.path.exists('output'): shutil.rmtree('output')
    os.mkdir('output')

    # Submit the benchmark to the cluster.

    for pdb in pdbs:
        if debug:
            command = (
                    'qsub', '-q', 'short.q', '-l', 'h_rt=0:30:00',
                    'benchmark.py', pdb, '--fast'
            )
        else:
            command = (
                    'qsub', '-t', '1-500',
                    'benchmark.py', pdb,
            )

        print ' '.join(command)
        subprocess.call(command)


if __name__ == '__main__':
    import optparse

    # Use optparse because it's available on chef.

    parser = optparse.OptionParser(usage='%prog [options] <benchmark>')
    parser.add_option('--compile-only', '-c', action='store_true', dest='compile_only')
    parser.add_option('--execute-only', '-x', action='store_true', dest='execute_only')
    options, arguments = parser.parse_args()
    benchmark = set(arguments)

    if not arguments:
        print 'Usage: launch.py [options] <full|test|mini|pdbs...>'
        print
        print 'launch.py: error: must specify a benchmark to launch.'
        sys.exit(1)

    # Compile rosetta.

    check_hostname()

    if not options.execute_only:
        error_code = compile_rosetta()
        if error_code != 0:
            sys.exit(error_code)

    if options.compile_only:
        sys.exit(1)

    # Run a test job, if requested.

    if 'test' in arguments:
        pdbs = ['structures/1srp.pdb']
        debug = True

    # Run the full benchmark on the specified structures, otherwise.

    else:
        pdbs = set()
        debug = False

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

    error_code = run_benchmark(pdbs, debug=debug)
    sys.exit(error_code)

