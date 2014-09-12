#!/usr/bin/env python2

#$ -S /usr/bin/python
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -cwd

import os
import sys; sys.path.append(os.getcwd())
import optparse
import subprocess
import re

from libraries import utilities
from libraries import settings; settings.load(interactive=False)
from libraries import database

# Parse arguments.

usage = 'SGE_TASK_ID=<id> loop_benchmark.py [options] <benchmark_id> <script>'
parser = optparse.OptionParser(usage=usage)
parser.add_option('--var', dest='vars', action='append', default=[])
parser.add_option('--fast', dest='fast', action='store_true')
options, arguments = parser.parse_args()

if len(arguments) != 2:
    print 'Usage:', usage
    print 
    print 'benchmark.py: error: expected 2 positional argument, got {0}.'.format(len(arguments))
    sys.exit(1)

task_id = int(os.environ['SGE_TASK_ID']) - 1
benchmark_id = int(arguments[0])
script_path = arguments[1]

# Figure out which loop to benchmark.

with database.connect() as session:
    benchmark = session.query(database.Benchmarks).get(benchmark_id)
    input_pdbs = benchmark.input_pdbs
    pdb_path = input_pdbs[task_id % len(input_pdbs)].pdb_path

loop_path = re.sub('\.pdb(\.gz)?$', '.loop', pdb_path)

# Set LD_LIBRARY_PATH so that the MySQL libraries can be found.

rosetta_env = os.environ.copy()
mysql_lib = '/netapp/home/kbarlow/lib/mysql-connector-c-6.1.2-linux-glibc2.5-x86_64/lib:'

try:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib + ':' + rosetta_env['LD_LIBRARY_PATH']
except KeyError:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib

# Build the RosettaScripts command line.

rosetta_path = os.path.abspath(settings.rosetta)
rosetta_scripts = os.path.join(rosetta_path, 'source', 'bin', 'rosetta_scripts.mysql.linuxgccrelease')
rosetta_database = os.path.join(rosetta_path, 'database')

rosetta_command = [
        rosetta_scripts,
        '-database', rosetta_database,
        '-in:file:s', pdb_path,
        '-in:file:native', pdb_path,
        '-inout:dbms:mode', 'mysql',
        '-inout:dbms:database_name', settings.db_name,
        '-inout:dbms:user', settings.db_user,
        '-inout:dbms:password', settings.db_password,
        '-inout:dbms:host', settings.db_host,
        '-inout:dbms:port', settings.db_port,
        '-out:nooutput',
        '-parser:protocol', script_path,
        '-parser:script_vars',
            'loop_file={0}'.format(loop_path),
            'fast={0}'.format('yes' if options.fast else 'no'),
]         + options.vars

# Run the benchmark.

process = subprocess.Popen(rosetta_command, env=rosetta_env,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

stdout, stderr = process.communicate(); process.poll()
sys.stdout.write(stdout); sys.stderr.write(stderr)
return_code = process.returncode

# Create a mapping between the benchmark and the protocol in the database.  

protocol_match = re.search("protocol_id '([1-9][0-9]*)'", stdout)
protocol_id = protocol_match.groups()[0] if protocol_match else None

with database.connect() as session:
    log_row = database.TracerLogs(benchmark_id, protocol_id, stdout, stderr)
    session.add(log_row)

    if protocol_id is not None:
        benchmark_map = database.BenchmarkProtocols(benchmark_id, protocol_id)
        session.add(benchmark_map)

sys.exit(return_code)
