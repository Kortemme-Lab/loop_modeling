#!/usr/bin/env python

import subprocess
from . import settings
from . import utilities

def install_dependencies_if_necessary():
    install_sqlalchemy()
    install_mysql_connector()
    settings.install()

def install_sqlalchemy():
    try:
        import sqlalchemy
    except ImportError:
        utilities.require_chef()
        raw_input("Installing sqlalchemy. Press [enter] to continue.")

        root_dir = utilties.get_benchmark_root()
        deps_dir = os.path.join(root, 'dependencies')
        package_dir = os.path.join(deps_dir, 'SQLAlchemy-0.9.7')
        package_archive = package_dir + '.tar.gz'

        unpack_command = 'tar', '-xz', '-C', deps_dir, '-f', package_archive
        install_command = '; '.join(
                'cd {0}'.format(parckage_dir),
                'python setup.py install --user',
        )
        install_command_32_bit = 'ssh', 'xeonint', install_command
        install_command_64_bit = 'ssh', 'iqint', install_command

        subprocess.check_call(unpack_command)
        subprocess.check_call(install_command_32_bit)
        subprocess.check_call(install_command_64_bit)

        print
    
def install_mysql_connector():
    try:
        import mysql.connector
    except ImportError:
        utilities.require_chef()
        raw_input("Installing mysql-connector. Press [enter] to continue.")

        root_dir = utilties.get_benchmark_root()
        deps_dir = os.path.join(root, 'dependencies')
        package_dir = os.path.join(deps_dir, 'mysql-connector-python-1.2.2')
        package_archive = package_dir + '.zip'

        unpack_command = 'unzip', '-d', deps_dir, package_archive
        install_command = shlex.split('python setup.py install --user')

        subprocess.check_call(unpack_command)
        subprocess.check_call(install_command, cwd=package_dir)

        print
