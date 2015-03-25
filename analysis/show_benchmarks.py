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
