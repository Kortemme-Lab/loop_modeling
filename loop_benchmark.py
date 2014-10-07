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
import json

from libraries import utilities
from libraries import settings; settings.load(interactive=False)
from libraries import database

# Parse arguments.

if len(sys.argv) != 2 or 'SGE_TASK_ID' not in os.environ:
    print 'Usage: SGE_TASK_ID=<id> loop_benchmark.py <benchmark_id>'
    sys.exit(1)

task_id = int(os.environ['SGE_TASK_ID']) - 1
benchmark_id = int(sys.argv[1])

# Figure out which loop to benchmark.

with database.connect() as session:
    benchmark = session.query(database.Benchmarks).get(benchmark_id)
    script_path = benchmark.rosetta_script
    script_vars = json.loads(benchmark.rosetta_script_vars or '[]')
    flags_path = benchmark.rosetta_flags
    fast = benchmark.fast
    input_pdbs = benchmark.input_pdbs
    pdb_path = input_pdbs[task_id % len(input_pdbs)].pdb_path
    loop_path = re.sub('\.pdb(\.gz)?$', '.loop', pdb_path)
    non_random = benchmark.non_random

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
            'fast={0}'.format('yes' if fast else 'no'),
]         + script_vars

if non_random:
    rosetta_command += ['-run:constant_seed', '-run:jran', task_id]

if flags_path is not None:
    rosetta_command += ['@', flags_path]

# Run the benchmark.

stdout, stderr = utilities.tee(rosetta_command, env=rosetta_env)

# Associate this run with the right benchmark and save log files.

protocol_match = re.search("protocol_id '([1-9][0-9]*)'", stdout)
protocol_id = protocol_match.groups()[0] if protocol_match else None

with database.connect() as session:
    if protocol_id is not None:
        benchmark_map = database.BenchmarkProtocols(benchmark_id, protocol_id)
        session.add(benchmark_map)
        session.commit()  # Make sure the protocol mapping is saved even if 
                          # something else messes up this transaction later on.
    log_row = database.TracerLogs(benchmark_id, protocol_id, stdout, stderr)
    session.add(log_row)
