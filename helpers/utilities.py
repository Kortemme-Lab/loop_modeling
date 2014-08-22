#!/usr/bin/env python

def get_benchmark_root():
    command = shlex.split('git rev-parse --show-toplevel')
    directory = subprocess.check_output(command)
    return directory.strip()

def require_chef():
    from socket import gethostname
    if gethostname() != 'chef.compbio.ucsf.edu':
        raise SystemExit("This script must be run on chef.")

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

