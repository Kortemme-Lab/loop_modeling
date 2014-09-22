#!/usr/bin/env python2

"""\
Display some identifying information for each benchmark in the database.

Usage: show_benchmarks.py
"""

import textwrap

from libraries import settings; settings.load()
from libraries import database
from libraries import docopt

docopt.docopt(__doc__)

with database.connect() as session:
    benchmarks = session.query(database.Benchmarks).all()

    columns = '#', 'Name', 'Description'
    column_widths = map(len, columns)

    for benchmark in benchmarks:
        column_widths = [
                max(column_widths[0], len(str(benchmark.id))),
                max(column_widths[1], len(benchmark.name or '???')),
                max(column_widths[2], len(benchmark.description or '')),
        ]

    total_width = sum(column_widths)
    column_widths[-1] -= max(total_width - 79 + 2, 0)

    row_template = '{{0:<{0[0]}}} {{1:<{0[1]}}} {{2}}'.format(
            column_widths)

    print row_template.format(*columns)
    print ' '.join('=' * x for x in column_widths)

    for benchmark in benchmarks:
        description = textwrap.fill(
                benchmark.description or '',
                width=column_widths[-1],
                initial_indent='',
                subsequent_indent=' ' * (79 - column_widths[-1])
        )
        print row_template.format(
                benchmark.id,
                benchmark.name or '???',
                description
        )
