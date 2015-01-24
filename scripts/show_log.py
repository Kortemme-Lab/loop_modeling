#!/usr/bin/env python2
# This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

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

