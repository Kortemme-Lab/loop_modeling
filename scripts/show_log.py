#!/usr/bin/env python2

"""\
Display the stdout and stderr logs for a particular benchmark run.

Usage:
    show_logs.py <benchmark_id> [<protocol_id>]
"""

from libraries import settings; settings.load()
from libraries import database
from libraries import docopt

arguments = docopt.docopt(__doc__)
benchmark_id = arguments['<benchmark_id>']
protocol_id = arguments['<protocol_id>']

with database.connect() as session:
    query = session.query(database.TracerLogs).\
            filter_by(benchmark_id=benchmark_id)

    if protocol_id is not None:
        query = query.filter_by(protocol_id=protocol_id)

    tracer_log = query.first()

    if tracer_log is None:
        print "No log found."
        raise SystemExit

    if tracer_log.stdout:
        print tracer_log.stdout
    if tracer_log.stderr:
        print tracer_log.stderr
