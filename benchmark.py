#!/usr/bin/env python

#$ -S /usr/bin/python
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=03:00:00
#$ -o /dev/null
#$ -j y
#$ -cwd

import os
import sys; sys.path.append(os.getcwd())
import optparse
import subprocess
import re

from helpers import utilities
from helpers import settings; settings.load(interactive=False)
from helpers import database

# Parse arguments (e.g. input pdb)

parser = optparse.OptionParser(usage='%prog [options] <pdb>')
parser.add_option('--id', dest='id', type=int, default=0)
parser.add_option('--var', dest='vars', action='append', default=[])
parser.add_option('--fast', dest='fast', action='store_true')
options, arguments = parser.parse_args()

if len(arguments) != 2:
    print 'Usage: benchmark.py [options] <script> <pdb>'
    print 
    print 'benchmark.py: error: expected 2 positional argument, got {0}.'.format(len(arguments))
    sys.exit(1)

script_path = arguments[0]
pdb_path = arguments[1]
loop_path = re.sub('\.pdb(\.gz)?$', '.loop', pdb_path)

# Setup LD_LIBRARY_PATH

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

process = subprocess.Popen(rosetta_command, env=rosetta_env, stdout=subprocess.PIPE)
stdout, stderr = process.communicate(); process.poll()
return_code = process.returncode

# Create a mapping between the benchmark and the protocol in the database.

protocol_match = re.search("protocol_id '([1-9][0-9]*)'", stdout)

if protocol_match:
    benchmark_id = options.id
    protocol_id = protocol_match.groups()[0]

    with database.connect() as session:
        benchmark_map = database.BenchmarkProtocols(benchmark_id, protocol_id)
        tracer_output = database.TracerOutput(protocol_id, stdout)
        session.add_all([benchmark_map, tracer_output])

sys.exit(return_code)
