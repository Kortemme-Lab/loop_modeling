#$ -S /usr/bin/python
#$ -l mem_free=2G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=03:00:00
#$ -cwd

# Setup LD_LIBRARY_PATH

rosetta_env = os.environ.copy()
mysql_lib = '/netapp/home/kbarlow/lib/mysql-connector-c-6.1.2-linux-glibc2.5-x86_64/lib:'

try:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib + ':' + rosetta_env['LD_LIBRARY_PATH']
except KeyError:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib

# Build the RosettaScripts command line.

rosetta_path = os.path.abspath(conf.path)
rosetta_scripts = os.path.join(rosetta_path, 'source', 'bin', 'rosetta_scripts')
rosetta_database = os.path.join(rosetta_path, 'database')

rosetta_command = (
        rosetta_scripts,
        '-database', rosetta_database,
        '-in:file:s', pdb_path,
        '-in:file:native', pdb_path,
        '-inout:database_mode', 'mysql',
        '-inout:database_filename', conf.db_name,
        '-mysql:host', conf.db_host,
        '-mysql:user', conf.db_user,
        '-mysql:password', conf.db_password,
        '-mysql:port', conf.db_port,
        '-out:use_database'
        '-parser:protocol', 'loop_modeler.xml',
        '-parser:script_vars', 'model={}'.format(pdb_tag),
)
# Run the benchmark.

error_code = subprocess.call(rosetta_command, env=rosetta_env)
sys.exit(error_code)
