#!/usr/bin/env python2

import os.path
import subprocess
import shlex
from . import settings
from . import utilities

def ask_to_install(message):
    try:
        raw_input("{}  Press [enter] to continue. ".format(message))
    except KeyboardInterrupt:
        print
        print "Aborting because required libary not installed."
        raise SystemExit


def require_sqlalchemy():
    try:
        import sqlalchemy

    except ImportError:
        utilities.require_chef()
        ask_to_install("Installing sqlalchemy.")

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
    
def require_mysql_connector():
    try:
        import mysql.connector

    except ImportError:
        utilities.require_chef()
        ask_to_install("Installing mysql-connector.")

        root_dir = utilties.get_benchmark_root()
        deps_dir = os.path.join(root, 'dependencies')
        package_dir = os.path.join(deps_dir, 'mysql-connector-python-1.2.2')
        package_archive = package_dir + '.zip'

        unpack_command = 'unzip', '-d', deps_dir, package_archive
        install_command = shlex.split('python setup.py install --user')

        subprocess.check_call(unpack_command)
        subprocess.check_call(install_command, cwd=package_dir)

        print

def require_pandas():
    try:
        message = "Installing pandas"
        import pandas
        major_version = int(pandas.__version__.split('.')[1])
        if major_version < 11:
            message = "Upgrading old version of pandas."
            raise ImportError

    except ImportError:
        ask_to_install(message)

        # Since this isn't meant to be used on chef, we can use pip to handle 
        # avoid to automatically pull down other dependencies if necessary.

        root_dir = utilities.get_benchmark_root()
        deps_dir = os.path.join(root_dir, 'dependencies')
        package_archive = os.path.join(deps_dir, 'pandas-0.14.1.zip')
        install_command = 'pip', 'install', '--user', package_archive

        subprocess.check_call(install_command)

        print

    import pandas
    return pandas

