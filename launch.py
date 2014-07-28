#!/usr/bin/env python2

import shutil
import subprocess
import conf

def check_hostname():
    from socket import gethostname

    if gethostname() != 'chef.compbio.ucsf.edu':
        raise SystemExit("This script must be run on chef.")

def compile_rosetta():
    rosetta_path = os.path.abspath(conf.path)

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
    
    for header in glob.glob(mysql_headers):
        shutil.copy(header, rosetta_headers)

    # Compile rosetta.

    scons_path = os.path.join(rosetta_path, 'source')

    compile_command = 'ssh', 'iqint', '; '.join([
            'cd "{}"'.format(scons_path),
            'nohup nice ./scons.py bin mode=release extras=mysql'])

    return subprocess.call(compile_command)

def test_benchmark():
    benchmark_path = os.getcwd()

    test_command = 'ssh', 'iqint', '; '.join([
        'cd "{}"'.format(benchmark_path),
        './benchmark.py structures/1srp.pdb --fast'])

    return subprocess.call(test_command)

def run_benchmark(pdbs):
    pdbs = [x for x in sorted(pdbs)]

    # Make sure all the inputs actually exist.

    for pdb in pdbs:
        if not os.path.exists(pdb):
            raise ValueError("'{}' does not exist.".format(pdb))

    # Submit the benchmark to the cluster.

    for pdb in pdbs:
        command = 'qsub', '-t', '1-500', 'benchmark.py', pdb
        subprocess.call(command)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('benchmark', nargs='+')
    parser.add_argument('--compile-only', '-c', action='store_true')
    arguments = parser.parse_args()

    # Compile rosetta.

    check_hostname()
    error_code = compile_rosetta()

    if arguments.compile_only or error_code != 0:
        sys.exit(error_code)

    # Run a test job, if requested.

    if 'test' in arguments.benchmark:
        error_code = test_benchmark()
        sys.exit(error_code)

    # Run the full benchmark on the specified structures, otherwise.

    else:
        pdbs = set()
        benchmark = set(arguments.benchmark)

        if 'full' in benchmark:
            benchmark.remove('full')
            pdbs += glob.glob('structures/*.pdb')

        if 'mini' in benchmark:
            benchmark.remove('mini')
            pdbs += ["structures/1c5e.pdb",
                     "structures/1cyo.pdb",
                     "structures/1exm.pdb",
                     "structures/2cpl.pdb",
                     "structures/3cla.pdb",
                     "structures/1oth.pdb"]

        pdbs += benchmarks
        error_code = run_benchmark(pdbs)
        sys.exit(error_code)

