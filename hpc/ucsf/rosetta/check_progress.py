#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Shane O'Connor
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
Checks the status of a loops benchmark run.

Usage:
    check_progress.py [<benchmark_name>] [options]

Arguments:
    benchmark_name
        The name for this benchmark.  It's ok for several benchmark runs to
        have the same name.  The analysis script will automatically pick the
        most recent one when given an ambiguous name.

Options:
    --database DATABASE -d DATABASE
        By default, the script will query the database defined in the settings.conf file. This setting allows other databases to be queried.

    --summary
        Omits the progress breakdown per structure and just reports the global progress. [default: 1]

    --list
        Lists all available benchmarks in the database (identified by name) but does not print any progress report.
"""

import sys
from libraries import settings
from libraries.dataController import DataController
from libraries import colortext


def exit(message):
    sys.exit(message + '\n')


#def get_benchmark_list_by_name(database_name):
#    with database.connect(db_name = database_name) as session:
#        return [r.name for r in session.execute('SELECT DISTINCT name from benchmarks ORDER BY benchmark_id DESC')]


def get_progress(data_controller, database_name, benchmark_name):
    return data_controller.get_progress(database_name, benchmark_name)


def get_progress_for_terminal(data_controller, database_name, benchmark_name, only_show_summary):
    '''Write progress to terminal.'''
    progress_data = get_progress(data_controller, database_name, benchmark_name)
    s = []
    if progress_data:
        progress, num_failed, nstruct = progress_data['Progress'], progress_data['FailureCount'], progress_data['nstruct']
        progress_fns = [colortext.mblue, colortext.mred, colortext.morange, colortext.myellow, colortext.mgreen]
        progress_fn = progress_fns[int(min(max(0, progress), 100)/25)]
        failed_cases_str = ''
        if num_failed:
            failed_cases_str = colortext.mred('\nFailed cases                : {0}'.format(num_failed))

        s = [progress_data['Messages'] or '']
        s.append('''Progress                                       : {0}
Number of cases                                : {1}{2}
Predictions per case                           : {3}
Total number of predictions                    : {4}
Completed predictions (extra jobs not counted) : {5}\n'''.format(progress_fn('%0.2f%%' % progress), progress_data['StructureCount'], failed_cases_str, nstruct, progress_data['TotalCount'], progress_data['CompletedCount']))

        if not only_show_summary:
            t = 'Completed jobs per case'
            s.append(t + '\n' + ('-' * len(t)))
            for k, v in sorted(progress_data['CountPerStructure'].iteritems()):
                progress_fn = progress_fns[int((float(min(max(0, v), nstruct))/float(nstruct))/0.25)]
                s.append('{0}: {1}'.format(k, progress_fn(v)))
            s.append('')
    return '\n'.join(s)


def report_progress(data_controller, database_name, benchmark_name, only_show_summary):
    print(get_progress_for_terminal(data_controller, database_name, benchmark_name, only_show_summary))


if __name__ == '__main__':
    try:
        from libraries import docopt
        arguments = docopt.docopt(__doc__)
        settings.load()
        benchmark_name = arguments['<benchmark_name>']
        summary = arguments['--summary']
        database_name = arguments['--database'] or settings.db_name
        data_controller = DataController(platform='database') if arguments['--database'] else DataController(platform='disk') 
        if arguments['--list']:
            name_list = data_controller.get_benchmark_list_by_name(database_name)
            if name_list:
                colortext.pgreen('\nList of available benchmarks in the {0} database:\n'.format(database_name))
                for n in name_list:
                    colortext.pcyan('  - {0}'.format(n))
                print('')
            else:
                colortext.pred('\nNo benchmarks defined in the {0} database.\n'.format(database_name))
        else:
            report_progress(data_controller, database_name, benchmark_name, summary)
    except KeyboardInterrupt:
        print


