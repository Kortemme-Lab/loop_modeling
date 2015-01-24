#!/usr/bin/env python2
# This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

import os.path
import subprocess
import shlex
from . import utilities

def ask_to_install(message):
    try:
        raw_input("{0}  Press [enter] to continue. ".format(message))
    except KeyboardInterrupt:
        print
        print "Aborting because required library not installed."
        raise SystemExit


def require_sqlalchemy():
    try:
        import sqlalchemy

    except ImportError:
        ask_to_install("Installing sqlalchemy.")

        root_dir = utilities.get_benchmark_root()
        libs_dir = os.path.join(root_dir, 'libraries')
        package_dir = os.path.join(libs_dir, 'SQLAlchemy-0.9.7')
        package_archive = package_dir + '.tar.gz'

        unpack_command = 'tar', '-xzv', '-C', libs_dir, '-f', package_archive
        install_command = '; '.join([
                'cd {0}'.format(package_dir),
                'python setup.py install --user',
        ])

        if utilities.is_this_chef():
            install_command_32_bit = 'ssh', 'xeonint', install_command
            install_command_64_bit = 'ssh', 'iqint', install_command

            subprocess.check_call(unpack_command)
            subprocess.check_call(install_command_32_bit)
            subprocess.check_call(install_command_64_bit)

        else:
            subprocess.check_call(unpack_command)
            subprocess.check_call(install_command, shell=True)

        print
    
def require_mysql_connector():
    try:
        import mysql.connector

    except ImportError:
        ask_to_install("Installing mysql-connector.")

        root_dir = utilities.get_benchmark_root()
        libs_dir = os.path.join(root_dir, 'libraries')
        package_dir = os.path.join(libs_dir, 'mysql-connector-python-1.2.2')
        package_archive = package_dir + '.zip'

        unpack_command = 'unzip', '-d', libs_dir, package_archive
        install_command = shlex.split('python setup.py install --user')

        subprocess.check_call(unpack_command)
        subprocess.check_call(install_command, cwd=package_dir)

        print

