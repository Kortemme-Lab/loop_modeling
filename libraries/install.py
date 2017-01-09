#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Kale Kundert
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

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
    except Exception:
        print('An exception occurred importing sqlalchemy.')
        raise
    
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


def require_biopython():
    try:
        import Bio.PDB
    except ImportError:
        ask_to_install("Installing biopython.")
        
        root_dir = utilities.get_benchmark_root()
        libs_dir = os.path.join(root_dir, 'libraries')
        package_dir = os.path.join(libs_dir, 'biopython-1.67')
        package_archive = package_dir + '.tar.gz'

        unpack_command = 'unzip', '-d', libs_dir, package_archive
        install_command = shlex.split('python setup.py install --user')

        subprocess.check_call(unpack_command)
        subprocess.check_call(install_command, cwd=package_dir)

        print


def require_flufl_lock():
    try:
        import flufl.lock
    except ImportError:
        ask_to_install("Installing flufl.lock")
        
        root_dir = utilities.get_benchmark_root()
        libs_dir = os.path.join(root_dir, 'libraries')
        package_dir = os.path.join(libs_dir, 'flufl.lock-2.0')
        package_archive = package_dir + '.tar.gz'

        unpack_command = 'tar', '-xzv', '-C', libs_dir, '-f', package_archive
        install_command = '; '.join([
                'cd {0}'.format(package_dir),
                'python2 setup.py install --user',
        ])

        subprocess.check_call(unpack_command)
        subprocess.check_call(install_command, shell=True)

        print


def require_klab():
    try:
        import klab
    except ImportError:
        ask_to_install("Installing klab.")

        install_command = ['pip', 'install', '--user', '-e',
                           'git+https://github.com/Kortemme-Lab/klab.git#egg=klab']
        subprocess.check_call(install_command)

        print
