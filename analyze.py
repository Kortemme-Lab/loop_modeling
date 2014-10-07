#!/usr/bin/python2

"""\
Create a report summarizing the given loops benchmark runs.  If multiple runs 
are given, extra plot explicitly comparing them will be included in the report.

Usage:
    analyze.py [options] <benchmarks>...

Options:
    -o FILE --output FILE
        File name for the PDF report being generated [default: report.pdf]

    --limit NUM -l NUM
        Limit the analysis to the given number of models per loop.  Usually 
        this is done to compare two or more benchmarks of different sizes.

    --keep-latex DIR
        Specify a directory where all the intermediate LaTeX and gnuplot 
        files should be saved.  These files can be used to integrate the report 
        into larger LaTeX documents.  By default these files are created in 
        /tmp and destroyed after the report is generated.

    {settings.config_args}

    --author AUTHOR
        Override the default author setting.  This name is used on the title 
        page of the generated report.

    {settings.database_args}

    --verbose
        Output progress messages and debugging information.

Authors:
    Roland A. Pache: Original implementation
    Kale Kundert: Database integration

Copyright:
    Copyright (C) 2011, 2012.  This software is not released under any license.
"""

# Imports (fold)
from __future__ import division

import collections
import contextlib
import math
import numpy
import os
import re
import shutil
import sys
import time

from libraries import colors
from libraries import settings
from libraries import statistics
from libraries import utilities

# Global variables (fold)
top_x = 5


class Report:

    @staticmethod
    def from_docopt_args(arguments):
        benchmarks = Benchmark.from_names(
                arguments['<benchmarks>'])

        report = Report(benchmarks)
        report.latex_dir = arguments['--keep-latex']
        report.keep_latex = arguments['--keep-latex'] is not None
        report.verbose = arguments['--verbose']

        # Throw out extra data if the user wants to limit the analysis to a 
        # certain number of structures.

        if arguments['--limit'] is not None:
            max_num_models = int(arguments['--limit'])

            for benchmark in benchmarks:
                for loop in benchmark:
                    loop.models = loop.models[:max_num_models]

        return report


    def __init__(self, benchmarks):
        self.benchmarks = benchmarks
        self.latex_dir = None
        self.keep_latex = False

        # If any benchmarks have duplicate names, de-duplicate them.

        name_counts = collections.Counter()
        name_postfixes = {}

        for benchmark in benchmarks:
            name_counts[benchmark.name] += 1

        for name in name_counts:
            name_count = name_counts[name]
            name_postfixes[name] = (
                    [''] if name_count == 1 else
                    ['_{}'.format(x + 1) for x in range(name_count)])

        for benchmark in benchmarks:
            benchmark.name += name_postfixes[benchmark.name].pop(0)

    def __len__(self):
        return len(self.benchmarks)

    def __iter__(self):
        return iter(self.benchmarks)

    def __reversed__(self):
        return reversed(iter(self.benchmarks))

    def make_report(self, path):
        with self.setup_latex_dir():
            self.setup_benchmark_colors()

            tex_path = os.path.join(self.latex_dir, 'report.tex')
            pdf_path = os.path.join(self.latex_dir, 'report.pdf')

            chapter_paths = []

            chapter_path = self.make_comparison_chapter()
            if chapter_path is not None:
                chapter_paths.append(chapter_path)

            for benchmark in self:
                chapter_path = self.make_benchmark_chapter(benchmark)
                chapter_paths.append(chapter_path)

            report_template = '''\
\\documentclass{{report}}

\\usepackage{{booktabs}}
\\usepackage{{graphicx}}
\\usepackage{{fullpage}}
\\usepackage{{hyperref}}

\\renewcommand{{\\textfraction}}{{0.00}}

\\hypersetup{{
    colorlinks,
    citecolor=black,
    filecolor=black,
    linkcolor=black,
    urlcolor=black
}}

\\begin{{document}}

\\title{{Loop Benchmark Report}}
\\author{{{author}}}
\\maketitle
\\tableofcontents

{chapters}

\\end{{document}}
'''
            author = settings.author
            chapters = '\n'.join(
                '\input{{{0}}}'.format(path)
                for path in chapter_paths)

            with open(tex_path, 'w') as file:
                file.write(report_template.format(**locals()))

            print "Running pdflatex..."

            if self.verbose:
                tex_command = 'pdflatex', 'report.tex'
            else:
                tex_command = 'pdflatex', '-interaction=batchmode', 'report.tex'

            utilities.run_command(tex_command, cwd=self.latex_dir, verbose=self.verbose)
            utilities.run_command(tex_command, cwd=self.latex_dir, verbose=self.verbose)
            shutil.copy(pdf_path, path)

    @contextlib.contextmanager
    def setup_latex_dir(self):
        """
        Create the directory where all the *.tex and *.pdf files generated for 
        this report will be stored.  This method is meant to be used in a 
        'with' statement, so that it can automatically clean up the directory 
        once the program terminates (either normally or via an exception).
        """

        import shutil
        import tempfile

        try:

            # Create the root latex directory.

            if self.latex_dir is not None:
                utilities.clear_directory(self.latex_dir)
            else:
                self.latex_dir = tempfile.mkdtemp('.loop_benchmark')

            self.latex_dir = os.path.abspath(self.latex_dir)

            # Create subdirectories for every benchmark and loop.

            for benchmark in self:
                benchmark.latex_dir = os.path.join(
                        self.latex_dir, 
                        benchmark.name + '_analysis',
                )
                os.mkdir(benchmark.latex_dir)

                for loop in benchmark:
                    loop.latex_dir = os.path.join(
                            benchmark.latex_dir, 
                            loop.pdb_id,
                    )
                    os.mkdir(loop.latex_dir)

            # Yield control to the rest of the program.

            yield

        finally:

            # Remove the latex directory if the user doesn't want to keep it.

            if not self.keep_latex and self.latex_dir is not None:
                shutil.rmtree(self.latex_dir)

    def setup_benchmark_colors(self):
        """
        Assign a color to each benchmark.  This creates a consistent color 
        theme throughout the entire report.  The colors are picked 
        from the Tango palette.  These colors are designed to look good 
        together, but will start to repeat if more than six benchmarks are 
        being compared.
        """

        for index, benchmark in enumerate(self):
            benchmark.color = colors.from_cycle(index)

    def make_comparison_chapter(self):
        """
        Create a latex file that organizes all the comparison plots into a 
        single chapter.  That latex file is the only output of this method.
        """

        if len(self) == 1:
            return None

        print "Making the comparison chapter..."

        # Make the comparison plots

        pdf_paths = self.make_comparison_plots()

        # Build the latex chapter up figure-by-figure.

        figure_template = '''\
\\section{{{0}}}
\\begin{{figure}}[h]
\\centering
\\includegraphics[angle=90,width=6in]{{{1}}}
\\end{{figure}}
\\pagebreak

'''
        tex_chapter = '''\
\\chapter{Benchmark Comparisons}

'''
        tex_chapter += figure_template.format(
                "Comparing percentage of subangstrom decoys",
                pdf_paths.percents_subangstrom)

        tex_chapter += figure_template.format(
                "Comparing best RMSD among lowest {top_x} scoring decoys".format(top_x=top_x),
                pdf_paths.best_top_x)

        tex_chapter += figure_template.format(
                "Comparing RMSD of lowest scoring decoys",
                pdf_paths.lowest_scores)

        tex_chapter += figure_template.format(
                "Comparing runtimes",
                pdf_paths.runtimes)

        # Save the latex chapter to disk.
        
        tex_path = os.path.join(self.latex_dir, 'comparison_plots.tex')

        with open(tex_path, 'w') as file:
            file.write(tex_chapter)

        return tex_path

    def make_comparison_plots(self):
        class Paths: pass
        pdf_paths = Paths()

        # Make the 'percent subangstrom' comparison plot.

        distributions = collections.OrderedDict([
            (benchmark, benchmark.percents_subangstrom)
            for benchmark in self
        ])
        path_template = 'percentage_subangstrom_models'
        custom_gnuplot_commands = '''\
set ylabel "Fraction of sub-\305 models [%]" rotate by -90
'''
        pdf_paths.percents_subangstrom = self.make_comparison_plot(
                distributions, path_template, custom_gnuplot_commands) 

        # Make the 'best of the top X' comparison plot.

        distributions = collections.OrderedDict([
            (benchmark, [x.rmsd for x in benchmark.best_top_x_models])
            for benchmark in self
        ])
        path_template = 'best_models_rmsd'
        custom_gnuplot_commands = '''\
set ylabel "r.m.s. deviation to crystal loop [\305]" rotate by -90
f(x)=1
'''
        custom_plot_arguments = 'f(x) with lines ls 9 notitle'
        pdf_paths.best_top_x = self.make_comparison_plot(
                distributions, path_template,
                custom_gnuplot_commands, custom_plot_arguments) 

        # Make the 'rmsd of lowest scoring decoys' comparison plot.

        distributions = collections.OrderedDict([
            (benchmark, [x.rmsd for x in benchmark.lowest_score_models])
            for benchmark in self
        ])
        path_template = 'lowest_score_models_rmsd'
        custom_gnuplot_commands = '''\
set ylabel "r.m.s. deviation to crystal loop [\305]" rotate by -90
f(x)=1
'''
        custom_plot_arguments = 'f(x) with lines ls 9 notitle'
        pdf_paths.lowest_scores = self.make_comparison_plot(
                distributions, path_template,
                custom_gnuplot_commands, custom_plot_arguments) 

        # Make the 'runtime' comparison plot.

        distributions = collections.OrderedDict([
            (benchmark, [x / 60 for x in benchmark.all_runtimes])
            for benchmark in self
        ])
        path_template = 'all_models_runtime'
        custom_gnuplot_commands = '''\
set ylabel "Runtime [min]" rotate by -90
'''
        pdf_paths.runtimes = self.make_comparison_plot(
                distributions, path_template, custom_gnuplot_commands) 

        return pdf_paths

    def make_comparison_plot(self, distributions, path_template,
            custom_gnuplot_commands, custom_plot_arguments=''):
        """
        Create a plot comparing the same distribution from several different 
        benchmarks.  Examples include the percent subangstrom distribution and 
        the RMSDs of the lowest scoring decoys.  The resulting plot will have a 
        nicely colored box plot for each benchmark.

        Inputs
        ------
        distributions:
            A dictionary mapping benchmark objects to some sort of 
            distribution.  The box plot will be created using the distribution 
            and labeled using the benchmark object.

        path_template:
            The base file name used to create the TSV, GNU, and EPS files 
            generated by this method.

        custom_gnuplot_commands:
            A string containing custom commands to pass to gnuplot immediately 
            before the 'plot' command.  This is meant to be used for doing 
            things like labeling the axes or adding useful vertical lines.

        Outputs
        -------
        This method creates three files: a TSV file containing the raw data 
        being plotted, a GNU file containing the gnuplot commands used to 
        generate the plot, and an EPS file containing the plot itself.  The 
        path to the generated EPS file is also returned.
        """

        tsv_path = os.path.join(self.latex_dir, path_template + '.tsv')
        gnu_path = os.path.join(self.latex_dir, path_template + '.gnu')
        pdf_path = os.path.join(self.latex_dir, path_template + '.pdf')

        # Write data to TSV file that can be easily parsed by gnuplot.

        boxplot_header = '#' + '\t'.join([
                'Protocol',
                'x',
                'lower',
                'first_quartile',
                'median',
                'third_quartile',
                'upper']) + '\n'

        boxplot_row = '\t'.join([
                '{benchmark.name}',
                '{gnuplot_index}',
                '{stats.lower_whisker}',
                '{stats.first_quartile}',
                '{stats.median}',
                '{stats.third_quartile}',
                '{stats.upper_whisker}']) + '\n'

        outlier_header = '#' + '\t'.join([
                'Protocol',
                'x',
                'outlier']) + '\n'

        outlier_row = '\t'.join([
                '{benchmark.name}',
                '{gnuplot_index}',
                '{outlier}']) + '\n'

        for index, benchmark in enumerate(reversed(distributions)):
            distribution = {1: distributions[benchmark]}
            boxplots = statistics.tukeyBoxAndWhisker(distribution)
            stats, outliers = boxplots[1]
            gnuplot_index = index + 1

            if not outliers:
                outliers = '?'

            with open(tsv_path, 'a') as file:
                file.write(boxplot_header)
                file.write(boxplot_row.format(**locals()))
                file.write('\n\n')
                file.write(outlier_header)
                for outlier in outliers:
                    file.write(outlier_row.format(**locals()))
                file.write('\n\n')

        # Generate plot using gnuplot.

        x_range = len(distributions) + 1
        x_ticks = ', '.join([
            '"{0.title}" {1}'.format(benchmark, i+1)
            for i, benchmark in enumerate(reversed(distributions))
        ])
        fig_height = min(1 + len(self), 5)

        gnuplot_script = '''\
set autoscale
set border 31
set tics out
set terminal pdf size {fig_height},6
set xtics ({x_ticks}) rotate by -90
set xtics nomirror
set ytics autofreq rotate by -90 center
set ytics nomirror
set noy2tics
set nox2tics

set style line 1 lt 1 lc rgb "dark-magenta" lw 2
set style line 2 lt 1 lc rgb "blue" lw 5 ps 1 pt 7
set style line 3 lt 1 lc rgb "forest-green" lw 2 ps 2 pt 13
set style line 4 lt 1 lc rgb "gold" lw 2 ps 1 pt 7
set style line 5 lt 1 lc rgb "red" lw 2 ps 2 pt 13
set style line 6 lt 1 lc rgb "black" lw 2
set style line 7 lt 1 lc rgb "dark-gray" lw 2
set style line 8 lt 1 lc rgb "gray" lw 2
set style line 9 lt 2 lc rgb "dark-gray" lw 5
set style fill solid 0.5

set boxwidth 0.75
set key below right
set xrange [0:{x_range}]
set encoding iso_8859_1
set notitle
unset xlabel
set yrange [0:]
set output "{pdf_path}"
{custom_gnuplot_commands}
plot {plot_arguments}
'''

        plot_template = ', \\\n     '.join([
                '"{tsv_path}" index {box_plot_index} using 2:4:3:7:6 with candlesticks whiskerbars lt 1 lc rgb "{color}" lw 5 notitle',
                '"{tsv_path}" index {box_plot_index} using 2:5:5:5:5 with candlesticks lt 1 lc rgb "black" lw 5 notitle',
                '"{tsv_path}" index {outliers_index} using 2:3 with points lt 1 lc rgb "{color}" lw 5 ps 0.5 pt 7 notitle',
        ])

        plot_arguments = ', \\\n     '.join([

                plot_template.format(
                    tsv_path=tsv_path,
                    box_plot_index=2*i,
                    outliers_index=2*i+1,
                    color=benchmark.color)

                for i, benchmark in enumerate(reversed(distributions))
        ])

        if custom_plot_arguments:
            plot_arguments += ', ' + custom_plot_arguments

        with open(gnu_path, 'w') as file:
            file.write(gnuplot_script.format(**locals()))

        utilities.run_gnuplot(gnu_path, verbose=self.verbose)

        return pdf_path

    def make_benchmark_chapter(self, benchmark):
        print "Making the {0.title} chapter...".format(benchmark)

        # Calculate the median percent subangstrom decoys.

        median_percent_subA = numpy.median(benchmark.percents_subangstrom)

        # Calculate the average number of models per loop.

        models_per_loop = numpy.mean([x.num_models for x in benchmark])

        # Calculate some basic statistics for different several different sets 
        # of models.
        
        summary_table_path = self.make_summary_stats_table(benchmark)

        # Make the 'RMSD of the best of the top 5 models' and 'fraction 
        # sub-angstrom' box plots.

        box_plot_paths = self.make_summary_box_plots(benchmark)
        rmsd_plot_path = box_plot_paths[0]
        subA_plot_path = box_plot_paths[2]

        # Make the table summarizing benchmark performance on each loop.

        loop_stats_path = self.make_loop_stats_table(benchmark)

        # Make the figure showing score vs rmsd plots and rmsd histograms for 
        # each loop in the benchmark.

        loop_plots_path = self.make_loop_decoy_plots(benchmark)

        # Create the complete individual benchmark chapter.

        chapter_template = '''\
\\chapter{{{benchmark.title}}}
\\section{{Overall benchmark performance}}

\\begin{{figure}}[h]
    \\centering
    \\parbox{{3in}}{{\\includegraphics[width=3in]{{{rmsd_plot_path}}}}}
    \\parbox{{3in}}{{\\includegraphics[width=3in]{{{subA_plot_path}}}}}
\\end{{figure}}

\\begin{{table}}[h]
    \\centering
    \\input{{{summary_table_path}}}
\\end{{table}}

\\begin{{center}}
    Median fraction of sub-\\AA{{}} models: {median_percent_subA:.2f}\\%\\\\
    Average number of models per loop: {models_per_loop:.2f}
\\end{{center}}

\\pagebreak
\\section{{Individual results per loop}}

\\input{{{loop_stats_path}}}

\\input{{{loop_plots_path}}}
'''
        tex_path = os.path.join(benchmark.latex_dir, 'individual_report.tex')

        with open(tex_path, 'w') as file:
            file.write(chapter_template.format(**locals()))

        return tex_path

    def make_summary_stats_table(self, benchmark):
        median = lambda x: numpy.median(list(x))

        best_top_x_models = benchmark.best_top_x_models
        best_top_x_stats = {
                'label': 'Best of top {0} models'.format(top_x),
                'median_rmsd': median(x.rmsd for x in best_top_x_models),
                'median_score': median(x.score for x in best_top_x_models),
                'median_runtime': median(x.runtime for x in best_top_x_models),
        }

        lowest_score_models = benchmark.lowest_score_models
        lowest_score_stats = {
                'label': 'Lowest scoring models',
                'median_rmsd': median(x.rmsd for x in lowest_score_models),
                'median_score': median(x.score for x in lowest_score_models),
                'median_runtime': median(x.runtime for x in lowest_score_models),
        }

        lowest_rmsd_models = benchmark.lowest_rmsd_models
        lowest_rmsd_stats = {
                'label': 'Lowest RMSD models',
                'median_rmsd': median(x.rmsd for x in lowest_rmsd_models),
                'median_score': median(x.score for x in lowest_rmsd_models),
                'median_runtime': median(x.runtime for x in lowest_rmsd_models),
        }

        all_models = benchmark.all_models
        all_stats = {
                'label': 'All models',
                'median_rmsd': median(x.rmsd for x in all_models),
                'median_score': median(x.score for x in all_models),
                'median_runtime': median(x.runtime for x in all_models),
        }

        median_fraction_subangstrom = median(benchmark.percents_subangstrom)

        # Create the latex page.

        table_template = '''\
\\begin{{tabular}}{{lrrr}}
\\toprule
                          & Median rmsd & Median score & Median time \\\\
Model selection           &       (\AA) &        (REU) &       (sec) \\\\
\\midrule
{0[label]:25s} & {0[median_rmsd]:11.2f} & {0[median_score]:12.2f} & {0[median_runtime]:11.0f} \\\\
{1[label]:25s} & {1[median_rmsd]:11.2f} & {1[median_score]:12.2f} & {1[median_runtime]:11.0f} \\\\
{2[label]:25s} & {2[median_rmsd]:11.2f} & {2[median_score]:12.2f} & {2[median_runtime]:11.0f} \\\\
{3[label]:25s} & {3[median_rmsd]:11.2f} & {3[median_score]:12.2f} & {3[median_runtime]:11.0f} \\\\
\\bottomrule
\\end{{tabular}}
'''
        tex_table = table_template.format(
                best_top_x_stats, lowest_score_stats,
                lowest_rmsd_stats, all_stats)

        tex_path = os.path.join(benchmark.latex_dir, 'benchmark_stats.tex')

        with open(tex_path, 'w') as file:
            file.write(tex_table)

        return tex_path

    def make_summary_box_plots(self, benchmark):
        tsv_path = os.path.join(benchmark.latex_dir, 'best_model_dists.tsv')
        gnu_path = os.path.join(benchmark.latex_dir, 'best_model_dists.gnu')
        pdf_path_rmsd = os.path.join(benchmark.latex_dir, 'best_model_dists_rmsd.pdf')
        pdf_path_score = os.path.join(benchmark.latex_dir, 'best_model_dists_score.pdf')
        pdf_path_subA = os.path.join(benchmark.latex_dir, 'best_model_dists_percent_subangstrom.pdf')

        # Calculate box plot parameters.

        best_top_x_models = benchmark.best_top_x_models

        distributions = {
                1: [x.rmsd for x in best_top_x_models],
                2: [x.score for x in best_top_x_models],
                3: benchmark.percents_subangstrom,
        }

        box_plots = statistics.tukeyBoxAndWhisker(distributions)

        # Write box plot data to a tab-separated value (TSV) file that can 
        # easily be parsed by gnuplot.

        with open(tsv_path, 'w') as file:
            for x in box_plots:
                box_params, outliers = box_plots[x]

                file.write('#x\t'+'lower\t'+'first_quartile\t'+'median\t'+'third_quartile\t'+'upper\n')
                for item in box_params: file.write('{0}\t'.format(item))
                file.write('\n\n\n')
                file.write('#x\toutlier\n')
                for outlier in outliers: file.write('{0}\t{1}\n'.format(x, outlier))
                if not outliers: file.write('{0}\t?\n'.format(x))
                file.write('\n\n')

        # Write the gnuplot script and generate the EPS plots.

        gnuplot_script = '''\
set autoscale
set border 31
set tics out
set terminal pdf
set size ratio 1
set noxtics
set xrange [0.5:1.5]
set nox2tics
set ytics 1
set ytics nomirror
set noy2tics

set style line 1 lt 1 lc rgb "dark-magenta" lw 2
set style line 2 lt 1 lc rgb "{benchmark.color}" lw 5 pt 7
set style line 3 lt 1 lc rgb "{benchmark.color}" lw 5
set style line 4 lt 1 lc rgb "gold" lw 2
set style line 5 lt 1 lc rgb "red" lw 5 pt 7
set style line 6 lt 1 lc rgb "black" lw 5
set style line 7 lt 1 lc rgb "dark-gray" lw 2
set style line 8 lt 1 lc rgb "gray" lw 2
set style line 9 lt 0 lc rgb "black" lw 5

set boxwidth 0.25
set key tmargin
set title "Best models performance distribution"
set noxlabel
set style fill solid 0.5
set encoding iso_8859_1
set ylabel "r.m.s. deviation to crystal loop [\305]"
set output "{pdf_path_rmsd}"
f(x)=1
plot "{tsv_path}" index 0 using 1:3:2:6:5 with candlesticks whiskerbars ls 2 notitle axes x1y1,\\
     "{tsv_path}" index 0 using 1:4:4:4:4 with candlesticks ls 6 notitle,\\
     "{tsv_path}" index 1 using 1:2 with points ls 2 ps 0.5 pt 7 notitle,\\
     f(x) with lines ls 9 notitle

set ylabel "Rosetta all-atom score"
set xrange [1.5:2.5]
set ytics autofreq
set output "{pdf_path_score}"
plot "{tsv_path}" index 2 using 1:3:2:6:5 with candlesticks whiskerbars ls 5 notitle axes x1y1,\\
     "{tsv_path}" index 2 using 1:4:4:4:4 with candlesticks ls 6 notitle,\\
     "{tsv_path}" index 3 using 1:2 with points ls 5 ps 0.5 pt 7 notitle

set title "Protocol performance distribution"
set ylabel "Fraction sub-\305 models [%]"
set xrange [2.5:3.5]
set ytics 10
set output "{pdf_path_subA}"
plot "{tsv_path}" index 4 using 1:3:2:6:5 with candlesticks whiskerbars ls 3 notitle axes x1y1,\\
     "{tsv_path}" index 4 using 1:4:4:4:4 with candlesticks ls 6 notitle,\\
     "{tsv_path}" index 5 using 1:2 with points ls 3 ps 0.5 pt 7 notitle
'''
        with open(gnu_path, 'w') as file:
            file.write(gnuplot_script.format(**locals()))

        # If there are no outliers, gnuplot will produce a warning.  This is a 
        # pretty common occurrence, and I think it's really bad to produce 
        # warning message for common occurrences.  So instead I opt to ignore 
        # stderr.  This is a little dangerous.  It would probably be better to 
        # suppress only the exact warning I know about.  But if we were really 
        # interested in doing things the right way, we would use matplotlib 
        # instead of gnuplot.

        with open(os.devnull) as devnull:
            utilities.run_gnuplot(gnu_path, stderr=devnull, verbose=self.verbose)

        return pdf_path_rmsd, pdf_path_score, pdf_path_subA

    def make_loop_stats_table(self, benchmark):
        """
        Create a table that uses a few different metrics to quantify how well 
        each loop in each benchmark was sampled.  The latex file produced by 
        this method contains nothing but this table and is meant to be included 
        in other documents (using \input).
        """

        # Define row and table templates that can easily be filled in.

        table_template = '''\
\\begin{{table}}[h]
    \\centering
    %\\footnotesize
    \\begin{{tabular}}{{lr|lrr|lrr}}
    \\toprule
    & &
        \multicolumn{{3}}{{l}}{{\emph{{Best of Top {top_x} Models}}}} &
        \multicolumn{{3}}{{|l}}{{\emph{{Closest Model}}}} \\\\
    Loop & Models &
        ID & RMSD & Score &
        ID & RMSD & Score \\\\
    \\midrule
{rows} \\\\
    \\bottomrule
    \\end{{tabular}}
\\end{{table}}
'''
        row_template = '    ' + ' & '.join([
            '{0.pdb_id}', '{0.num_models:3}',
            '{1.id:3}', '{1.rmsd:4.2f}', '{1.score:7.2f}',
            '{2.id:3}', '{2.rmsd:4.2f}', '{2.score:7.2f}',
        ])
        empty_row_template = '    ' + ' & '.join([
            '{0.pdb_id}', '{0.num_models:3}',
            '---', ' ---', '    ---',
            '---', ' ---', '    ---',
        ])

        # Create a row for each loop in the benchmark.

        table_rows = []

        for loop in sorted(benchmark, key=lambda x: x.pdb_id):
            if loop.has_data:
                row = row_template.format(
                    loop,
                    loop.best_top_x_model,
                    loop.lowest_rmsd_model
                )
            else:
                row = empty_row_template.format(loop)

            table_rows.append(row)

        tex_table = table_template.format(
                top_x=top_x, rows=' \\\\\n'.join(table_rows))

        # Write the output latex file.

        tex_path = os.path.join(benchmark.latex_dir, 'loop_stats.tex')

        with open(tex_path, 'w') as file:
            file.write(tex_table)

        return tex_path

    def make_loop_decoy_plots(self, benchmark):
        pages = []
        page_template = '''\
\\begin{{figure}}[h]
\\centering
\\begin{{tabular}}{{cc}}
{0} \\\\
\\end{{tabular}}\n\
\\end{{figure}}
'''
        rows = []
        rows_per_page = 5
        row_template = (
                '\\includegraphics[height=1.7in]{{{0}}} & '
                '\\includegraphics[height=1.7in]{{{1}}}')
        missing_data_template = (
                '\\multicolumn{{2}}{{c}}{{'
                '\\parbox[t][1.7in][c]{{4in}}{{'
                '\\centering{{}}'
                'No data for {0.pdb_id}}}}}')

        for i, loop in enumerate(sorted(benchmark, key=lambda x: x.pdb_id)):
            score_vs_rmsd_paths = self.make_score_vs_rmsd_plot(loop)
            rmsd_histogram_path = self.make_rmsd_histogram(loop)

            if loop.has_data:
                row = row_template.format(
                        score_vs_rmsd_paths[1], rmsd_histogram_path)
            else:
                row = missing_data_template.format(loop)

            rows.append(row)

            if i % rows_per_page == rows_per_page - 1:
                page = page_template.format(' \\\\\n'.join(rows))
                pages.append(page)
                rows = []

        page = page_template.format(' \\\\\n'.join(rows))
        pages.append(page)

        tex_path = os.path.join(benchmark.latex_dir, 'loop_plots.tex')

        with open(tex_path, 'w') as file:
            file.write('\n'.join(pages))

        return tex_path

    def make_score_vs_rmsd_plot(self, loop):
        """
        Create a score vs RMSD plot for the given loop.  In fact two plots are 
        made: one which includes every model and one which includes only the 
        top 75% best scoring models.  Normally the second plot is of more 
        interest, because it focuses better on the interesting lower-left 
        region of the plot.  The full plots often have outliers that really 
        scale the score axis.
        """

        # This method would be much more concise if it used matplotlib.

        if not loop.has_data:
            return

        tsv_path = os.path.join(loop.latex_dir, 'score_vs_rmsd.tsv')
        gnu_path = os.path.join(loop.latex_dir, 'score_vs_rmsd.gnu')
        pdf_path_100 = os.path.join(loop.latex_dir, 'score_vs_rmsd_all.pdf')
        pdf_path_75 = os.path.join(loop.latex_dir, 'score_vs_rmsd_third_quartile.pdf')

        tsv_row = '{0.id}\t{0.rmsd}\t{0.score}\n'

        sorted_models = loop.models_sorted_by_score
        scores = loop.scores
        min_score, max_score = min(scores), max(scores)
        third_quartile = numpy.percentile(scores, 75)
        native_score = 0    # This isn't stored in the database yet.

        # Write score vs RMSD data to a tab-separated value (TSV) file that can 
        # easily be parsed by gnuplot.

        with open(tsv_path, 'w') as file:
            file.write('#Model\tLoop_rmsd\tTotal_score\n')
            file.write('input_structure\t0.0\t{0}\n'.format(native_score))

            # All models
            file.write('\n\n')
            for model in sorted_models:
                file.write(tsv_row.format(model))

            # Top X scoring models
            file.write('\n\n')
            for model in sorted_models[:top_x]:
                file.write(tsv_row.format(model))

            # Top scoring model
            file.write('\n\n')
            file.write(tsv_row.format(sorted_models[0]))

        # Write the gnuplot script and generate the EPS plots.

        gnuplot_script = '''\
set autoscale
set border 31
set tics out
set terminal pdf
set xtics autofreq
set xtics nomirror
set ytics autofreq
set ytics nomirror
set noy2tics
set nox2tics

set style line 1 lt 1 lc rgb "dark-magenta" lw 2
set style line 2 lt 1 lc rgb "{loop.benchmark.color}" lw 2 ps 0.5 pt 7
set style line 3 lt 1 lc rgb "forest-green" lw 2 ps 2 pt 13
set style line 4 lt 1 lc rgb "dark-gray" lw 2 ps 0.5 pt 7
set style line 5 lt 1 lc rgb "black" lw 2 ps 0.8 pt 13
set style line 6 lt 1 lc rgb "black" lw 2
set style line 7 lt 1 lc rgb "dark-gray" lw 2
set style line 8 lt 1 lc rgb "gray" lw 2
set style line 9 lt 2 lc rgb "dark-gray" lw 5

set boxwidth 0.75

set key below right
set xrange [0:]
set encoding iso_8859_1
set title "{loop.pdb_id}: {loop.percent_subangstrom:0.2f}% sub-\305 models"
set xlabel "r.m.s. deviation to crystal loop [\305]"
set arrow from 1, graph 0 to 1, graph 1 ls 9 nohead
set ylabel "Rosetta all-atom score"
set output "{pdf_path_100}"
plot "{tsv_path}" index 1 using ($2):($3) with points ls 2 title "all models" axes x1y1, \\
     "{tsv_path}" index 2 using ($2):($3) with points ls 4 title "5 lowest energy models" axes x1y1, \\
     "{tsv_path}" index 3 using ($2):($3) with points ls 5 title "top 5 best model" axes x1y1
set yrange [:{third_quartile}]
set output "{pdf_path_75}"
set xrange [0:]
plot "{tsv_path}" index 1 using ($2):($3) with points ls 2 title "75% lowest-scoring models" axes x1y1, \\
     "{tsv_path}" index 2 using ($2):($3) with points ls 4 title "5 lowest energy models" axes x1y1, \\
     "{tsv_path}" index 3 using ($2):($3) with points ls 5 title "top 5 best model" axes x1y1
'''
        with open(gnu_path, 'w') as file:
            file.write(gnuplot_script.format(**locals()))

        utilities.run_gnuplot(gnu_path, verbose=self.verbose)

        return pdf_path_100, pdf_path_75

    def make_rmsd_histogram(self, loop):
        """
        Create a smoothed RMSD histogram for the given loop.  100 bins are used 
        when making the plot, and the smoothing is done by gnuplot.
        """

        # This method would be much more concise if it used matplotlib.

        if not loop.has_data:
            return

        tsv_path = os.path.join(loop.latex_dir, 'rmsd_histogram.tsv')
        gnu_path = os.path.join(loop.latex_dir, 'rmsd_histogram.gnu')
        pdf_path = os.path.join(loop.latex_dir, 'rmsd_histogram.pdf')

        # Write histogram data to a tab-separated value (TSV) file that can 
        # easily be parsed by gnuplot.

        num_bins = 100
        histogram = statistics.histogram(loop.rmsds, num_bins)

        with open(tsv_path, 'w') as file:
            file.write('#All models\n')
            file.write('#RMSD\tFrequency\n')

            for rmsd, count in histogram:
                count = num_bins * count / len(loop)
                file.write('{0}\t{1}\n'.format(rmsd, count))

        # Write the gnuplot script and generate the EPS plot.

        gnuplot_script='''\
set autoscale
set border 31
set tics out
set terminal pdf enhanced color
set xtics autofreq
set xtics nomirror
set ytics autofreq
set ytics nomirror
set noy2tics
set nox2tics

set style line 1 lt 1 lc rgb "dark-magenta" lw 2
set style line 2 lt 1 lc rgb "{loop.benchmark.color}" lw 8 ps 1 pt 7
set style line 3 lt 1 lc rgb "forest-green" lw 2 ps 2 pt 13
set style line 4 lt 1 lc rgb "gold" lw 2 ps 1 pt 7
set style line 5 lt 1 lc rgb "red" lw 2 ps 2 pt 13
set style line 6 lt 1 lc rgb "black" lw 2
set style line 7 lt 1 lc rgb "dark-gray" lw 2
set style line 8 lt 1 lc rgb "gray" lw 2
set style line 9 lt 2 lc rgb "dark-gray" lw 5

set boxwidth 0.75

set key below right
set xrange [0:]
set encoding iso_8859_1
set title "{loop.pdb_id}: {loop.percent_subangstrom:0.2f}% sub-\305 models"
set xlabel "r.m.s. deviation to crystal loop [\305]"
set yrange [0:]
set arrow from 1, graph 0 to 1, graph 1 ls 9 nohead
set ylabel "Fraction of models [%]"
set output "{pdf_path}"
plot "{tsv_path}" index 0 using ($1):($2) smooth bezier with lines ls 2 title "all models" axes x1y1
'''
        with open(gnu_path, 'w') as file:
            file.write(gnuplot_script.format(**locals()))

        utilities.run_gnuplot(gnu_path, verbose=self.verbose)

        return pdf_path


class Benchmark:

    @staticmethod
    def from_names(names):
        benchmarks = []

        # The list of names may include several different kinds of name:
        #
        # 1. Benchmark id (primary key into database defined in settings file)
        # 2. Benchmark name (most recent benchmark run with given name in database)
        # 3. Path to benchmark data (stored in tabular flat file format)

        for name in names:
            if name.endswith('.results'):
                benchmark = Benchmark.from_flat_file(name)
            else:
                benchmark = Benchmark.from_database(name)

            if benchmark:
                benchmarks.append(benchmark)
            else:
                message = "Skipping the empty {0} benchmark..."
                utilities.print_warning(message, benchmark.title)

        return benchmarks

    @staticmethod
    def from_flat_file(path, name=None):
        from os.path import basename, splitext

        name = name or splitext(basename(path))[0]
        benchmark = Benchmark(name)

        print "Loading the {0.title} benchmark from a flat file...".format(benchmark)

        with open(path) as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#'): continue

                tag, id, rmsd, score, runtime = line.split()
                id, rmsd, score, runtime = \
                        int(id), float(rmsd), float(score), int(runtime)

                loop = benchmark.loops.setdefault(tag, Loop(benchmark, tag))
                model = Model(loop, id, score, rmsd, runtime)
                loop.models.append(model)

        return benchmark

    @staticmethod
    def from_database(name_or_id):
        from libraries import database
        from libraries import settings

        with database.connect() as session:

            # Decide whether a name or id was used to specify a benchmark run, 
            # and load the corresponding data out of the database.  The meaning 
            # of name_or_id is inferred from its type: names are expected to be 
            # strings and ids are expected to be integers.  If more than one
            # benchmark has the same name, the most recent one will be used.

            try:
                id = int(name_or_id)
                db_benchmark = session.query(database.Benchmarks).get(id)

                if db_benchmark is None:
                    message = "No benchmark '{}' in the database."
                    utilities.print_error_and_die(message, id)

            except ValueError:
                name = name_or_id
                query = session.query(database.Benchmarks).filter_by(name=name)
                db_benchmark = sorted(query.all(), key=lambda x: x.start_time)[0]

            # Fill in the benchmark data structure from the database.
            
            benchmark = Benchmark(db_benchmark.name, db_benchmark.title)

            print "Loading the {0.title} benchmark from the database...".format(benchmark)

            for db_input in db_benchmark.input_pdbs:
                path = db_input.pdb_path
                benchmark.loops[path] = Loop(benchmark, path)

            for structure in db_benchmark.structures:
                loop = benchmark.loops[structure.input_tag]
                id = len(loop.models) + 1
                score = structure.score_features.score
                rmsd = structure.rmsd_features.protein_backbone
                runtime = structure.runtime_features.elapsed_time

                model = Model(loop, id, score, rmsd, runtime)
                loop.models.append(model)

        return benchmark


    def __init__(self, name, title=None):
        self.name = name
        self.manual_title = title
        self.loops = {}         # Set by Report.from_...()
        self.latex_dir = None   # Set by Report.setup_latex_dir()
        self.color = None       # Set by Report.setup_benchmark_colors()

    def __str__(self):
        return '<Benchmark name={0.name}>'.format(self)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __iter__(self):
        return self.loops.itervalues()

    def __len__(self):
        return len(self.loops)

    def __nonzero__(self):
        return any(self)

    def __getitem__(self, tag):
        return self.loops[tag]

    def __setitem__(self, tag, loop):
        self.loops[tag] = loop

    @property
    def title(self):
        if self.manual_title:
            return self.manual_title

        words = self.name.split('_')
        phrase = ' '.join(words)
        phrase = phrase.title()
        phrase = re.sub(r'\b[Kk][Ii][Cc]\b', 'KIC', phrase)
        phrase = re.sub(r'\b[Cc][Cc][Dd]\b', 'CCD', phrase)
        phrase = re.sub(r'\b[Nn][Gg][Kk]\b', 'NGK', phrase)
        return phrase

    @property
    def all_models(self):
        return sum((loop.models for loop in self if loop.has_data), [])

    @property
    def all_runtimes(self):
        return sum((loop.runtimes for loop in self if loop.has_data), [])

    @property
    def best_top_x_models(self):
        return [loop.best_top_x_model for loop in self if loop.has_data]

    @property
    def lowest_score_models(self):
        return [loop.lowest_score_model for loop in self if loop.has_data]

    @property
    def lowest_rmsd_models(self):
        return [loop.lowest_rmsd_model for loop in self if loop.has_data]

    @property
    def percents_subangstrom(self):
        return [loop.percent_subangstrom for loop in self if loop.has_data]


class Loop:

    def __init__(self, benchmark, path):
        self.benchmark = benchmark
        self.path = path
        self.models = []        # Set by Report.from_...()
        self.latex_dir = None   # Set by Report.setup_latex_dir()

    def __iter__(self):
        return iter(self.models)

    def __len__(self):
        return len(self.models)

    def __nonzero__(self):
        return bool(self.models)

    @property
    def pdb_id(self):
        return os.path.basename(self.path)[0:4]

    @property
    def num_models(self):
        return len(self)

    @property
    def has_data(self):
        return len(self) > 0

    @property
    def scores(self):
        return [x.score for x in self.models]

    @property
    def rmsds(self):
        return [x.rmsd for x in self.models]

    @property
    def runtimes(self):
        return [x.runtime for x in self.models]

    @property
    def models_sorted_by_score(self):
        return sorted(self.models, key=lambda x: x.score)

    @property
    def models_sorted_by_rmsd(self):
        return sorted(self.models, key=lambda x: x.rmsd)

    @property
    def best_top_x_model(self):
        top_x_models = self.models_sorted_by_score[:top_x]
        return sorted(top_x_models, key=lambda x: x.rmsd)[0]

    @property
    def lowest_score_model(self):
        return self.models_sorted_by_score[0]

    @property
    def lowest_rmsd_model(self):
        return self.models_sorted_by_rmsd[0]

    @property
    def percent_subangstrom(self):
        num_sub_a = sum(1 for x in self.models if x.rmsd < 1.0)
        return 100 * num_sub_a / len(self.models)


class Model:

    def __init__(self, loop, id, score, rmsd, runtime):
        self.loop = loop
        self.id = id
        self.score = score
        self.rmsd = rmsd
        self.runtime = runtime



if __name__ == '__main__':
    from libraries import docopt

    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
        settings.load(arguments)
        report = Report.from_docopt_args(arguments)
        report.make_report(arguments['--output'])

    except KeyboardInterrupt:
        print

