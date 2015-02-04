#!/usr/bin/env python2
# This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

"""\
Display some identifying information for each benchmark in the database.

Usage:
    show_benchmarks.py [options]

Options:
    {settings.config_args}

    {settings.database_args}
"""

import textwrap

from libraries import settings
from libraries import database
from libraries import docopt

arguments = docopt.docopt(__doc__.format(**locals()))
settings.load(arguments)

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
    column_widths[-1] -= max(total_width - 79 + len(columns) - 1, 0)

    row_template = '{{0:<{0[0]}}} {{1:<{0[1]}}} {{2}}'.format(
            column_widths)

    print row_template.format(*columns)
    print ' '.join('=' * x for x in column_widths)

    for benchmark in benchmarks:
        indent = ' ' * (79 - column_widths[-1])
        description = textwrap.fill(
                benchmark.description or '',
                width=79,
                initial_indent=indent,
                subsequent_indent=indent,
        ).strip()

        print row_template.format(
                benchmark.id,
                benchmark.name or '???',
                description
        )
