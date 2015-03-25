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

"""\
Display the stdout and stderr logs for a particular benchmark run.

Usage:
    show_logs.py <benchmark_id> [options]

Options:
    --broken
        Filter out jobs that successfully produced a score vs rmsd point.

    --pdb PDB
        Filter out jobs that don't match the given PDB.

    {settings.config_args}

    {settings.database_args}
"""

# Add options to pick a particular structure.
# Add options to pick given structure (i.e 1-500 index) or lowest scoring.
# Add options to show logs for jobs that died.
# Default: only benchmark id required.  Default structure is the lowest energy 
# model for the first PDB tag sorted alphabetically.

from libraries import settings
from libraries import database
from libraries import docopt

arguments = docopt.docopt(__doc__.format(**locals()))
benchmark_id = arguments['<benchmark_id>']

settings.load(arguments)

with database.connect() as session:
    query = session.query(database.TracerLogs).\
            filter_by(benchmark_id=benchmark_id)

    if arguments['--pdb']:
        print "This filter is not yet supported."
        raise SystemExit

    if arguments['--broken']:
        query = query.filter_by(protocol_id=0)

    tracer_log = query.first()

    if tracer_log is None:
        print "No log found."
        raise SystemExit

    if tracer_log.stdout:
        print tracer_log.stdout
    if tracer_log.stderr:
        print tracer_log.stderr

