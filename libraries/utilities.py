#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Roland Pache
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

def get_benchmark_root():
    import shlex
    import subprocess
    command = shlex.split('git rev-parse --show-toplevel')
    directory = subprocess.check_output(command)
    return directory.strip()

def clear_directory(directory):
    import os, shutil
    if os.path.exists(directory): shutil.rmtree(directory)
    os.mkdir(directory)
def is_this_chef():
    from socket import gethostname
    return gethostname() == 'chef.compbio.ucsf.edu'

def require_chef():
    if not is_this_chef():
        raise SystemExit("This script must be run on chef.")

def tee(*popenargs, **kwargs):
    import subprocess, select, sys

    process = subprocess.Popen(
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            *popenargs, **kwargs)

    stdout, stderr = '', ''

    print "Process ID:", process.pid

    def read_stream(input_callback, output_stream):   # (no fold)
        read = input_callback()
        output_stream.write(read)
        output_stream.flush()
        return read

    while process.poll() is None:
        watch = process.stdout.fileno(), process.stderr.fileno()
        ready = select.select(watch, [], [])[0]

        for fd in ready:
            if fd == process.stdout.fileno():
                stdout += read_stream(process.stdout.readline, sys.stdout)
            if fd == process.stderr.fileno():
                stderr += read_stream(process.stderr.readline, sys.stderr)

    stdout += read_stream(process.stdout.read, sys.stdout)
    stderr += read_stream(process.stderr.read, sys.stderr)

    return stdout, stderr

def check_output(*popenargs, **kwargs):
    """ 
    Run command with arguments and return its output as a byte string.

    If the exit code was non-zero it raises a CalledProcessError.  The
    CalledProcessError object will have the return code in the returncode
    attribute and output in the output attribute.

    The arguments are the same as for the Popen constructor.  Example:

    >>> check_output(["ls", "-l", "/dev/null"])
    'crw-rw-rw- 1 root root 1, 3 Oct 18  2007 /dev/null\n'

    The stdout argument is not allowed as it is used internally.
    To capture standard error in the result, use stderr=STDOUT.

    >>> check_output(["/bin/sh", "-c",
    ...               "ls -l non_existent_file ; exit 0"],
    ...              stderr=STDOUT)
    'ls: non_existent_file: No such file or directory\n'
    """
    from subprocess import Popen, PIPE, CalledProcessError
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = Popen(stdout=PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        raise CalledProcessError(retcode, cmd)
    return output

# Monkey patch the subprocess module to add the useful check_output() function.

import subprocess
subprocess.check_output = check_output

def run_command(cmd, **kwargs):
    verbose = kwargs.pop('verbose', False)
    cwd = kwargs.get('cwd', None)

    if verbose:
        if cwd is not None: print '$ cd {0}'.format(cwd)
        print '$', ' '.join(cmd)
        subprocess.check_call(cmd, **kwargs)
    else:
        subprocess.check_output(cmd, **kwargs)

def run_gnuplot(gnuplot_script, **kwargs):
    run_command(('gnuplot', gnuplot_script), **kwargs)

def print_warning(message, *args, **kwargs):
    if args or kwargs: message = message.format(*args, **kwargs)
    print '\033[1;31m' + message + '\033[0;0m'

def print_error_and_die(message, *args, **kwargs):
    print_warning(message + "  Aborting...", *args, **kwargs)
    raise SystemExit(1)

