#!/usr/bin/env python

#$ -S /usr/bin/python
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=03:00:00
#$ -o output
#$ -j y
#$ -cwd

import os
import sys; sys.path.append(os.getcwd())
import optparse
import subprocess
import re
import conf

# Parse arguments (e.g. input pdb)

parser = optparse.OptionParser(usage='%prog [options] <pdb>')
parser.add_option('--fast', action='store_true', dest='fast')
options, arguments = parser.parse_args()

if len(arguments) != 1:
    print 'Usage: benchmark.py [options] <pdb>'
    print 
    print 'benchmark.py: error: must specify exactly one input pdb.'
    sys.exit(1)

pdb_path = arguments[0]
loop_path = re.sub('\.pdb(\.gz)?$', '.loop', pdb_path)

# Setup LD_LIBRARY_PATH

rosetta_env = os.environ.copy()
mysql_lib = '/netapp/home/kbarlow/lib/mysql-connector-c-6.1.2-linux-glibc2.5-x86_64/lib:'

try:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib + ':' + rosetta_env['LD_LIBRARY_PATH']
except KeyError:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib

# Build the RosettaScripts command line.

rosetta_path = os.path.abspath(conf.rosetta)
rosetta_scripts = os.path.join(rosetta_path, 'source', 'bin', 'rosetta_scripts.mysql.linuxgccrelease')
rosetta_database = os.path.join(rosetta_path, 'database')

rosetta_command = (
        rosetta_scripts,
        '-database', rosetta_database,
        '-in:file:s', pdb_path,
        '-in:file:native', pdb_path,
        '-inout:dbms:mode', 'mysql',
        '-inout:dbms:database_name', conf.db_name,
        '-inout:dbms:user', conf.db_user,
        '-inout:dbms:password', conf.db_password,
        '-inout:dbms:host', conf.db_host,
        '-inout:dbms:port', conf.db_port,
        '-out:use_database',
        '-parser:protocol', 'loop_modeler.xml',
        '-parser:script_vars',
            'loop_file={0}'.format(loop_path),
            'fast={0}'.format('yes' if options.fast else 'no'),
)
# Run the benchmark.

error_code = subprocess.call(rosetta_command, env=rosetta_env)
sys.exit(error_code)
