#!/usr/bin/env python2
# encoding: utf-8

""" 
Generate a report comparing any number of benchmark runs.  The report will 
include the following figures for each benchmark being compared:

1. The distribution of the RMSD of the lowest scoring decoys.
2. The distribution of the fraction of sub-angstrom decoys sampled.
3. Score vs. RMSD plots for each structure.

Usage:
    analysis.py <benchmark_names_or_ids>...
"""

import matplotlib.pyplot as plt

from helpers import settings; settings.load()
from helpers import database
from helpers import colors
from helpers import install

pd = pandas = install.require_pandas()

def extract_benchmark_data():
    load_table = database.data_frame_from_table

    with database.connect() as session:
        benchmarks = load_table(session, database.Benchmarks)
        benchmark_protocols = load_table(session, database.BenchmarkProtocols)
        protocols = load_table(session, database.Protocols)
        batches = load_table(session, database.Batches)
        structures = load_table(session, database.Structures)
        rmsds = load_table(session, database.ProteinRmsdNoSuperposition)
        scores = load_table(session, database.TotalScores)

    scores_vs_rmsds = pd.merge(scores, rmsds, on='struct_id')
    scores_vs_rmsds = pd.merge(scores_vs_rmsds, structures, on='struct_id')
    scores_vs_rmsds = pd.merge(scores_vs_rmsds, batches, on='batch_id')
    scores_vs_rmsds = pd.merge(scores_vs_rmsds, benchmark_protocols, on='protocol_id')

    scores_vs_rmsds = scores_vs_rmsds[[
        'benchmark_id', 'protocol_id', 'input_tag',
        'score', 'protein_backbone']]
    scores_vs_rmsds.rename(columns={
        'protein_backbone': 'rmsd',
        'input_tag': 'pdb'}, inplace=True)

    return scores_vs_rmsds

def tag_from_pdb(pdb):
    import os.path
    tag = os.path.basename(pdb)
    if tag.endswith('.gz'): tag = tag[:-3]
    if tag.endswith('.pdb'): tag = tag[:-4]
    return tag


def plot_subangstrom_decoy_distributions():
    pass

def plot_lowest_score_rmsd_distributions():
    pass

def plot_score_vs_rmsd_funnels(benchmark_ids, benchmarks):
    # Confusing variable names: benchmarks doesn't mean the same thing here as 
    # it does in main.

    rows, cols = 4, 3
    figure, axes = plt.subplots(rows, cols, sharex='col')

    # Can't get labels centered...
    #plt.set_xlabel(u'Backbone RMSD (Å)')
    #plt.set_ylabel('Score (REU)')

    # Breaks if benchmarks have different inputs...
    for id in benchmark_ids:
        benchmark = benchmarks[benchmarks.benchmark_id == id]
        inputs = 1 * list(set(benchmark.pdb))

        for index, pdb in enumerate(sorted(inputs)):
            tag = tag_from_pdb(pdb)
            row, col = index // cols, index % cols
            axis = axes[row, col]

            print row, col, tag

            data = benchmark[benchmark.pdb == pdb]
            print data.head()
            data = data.sort('score')

            min_y = data.score.iloc[0]
            max_y = data.score.iloc[int(0.75 * len(data.rmsd))]
            min_y -= 0.1 * (max_y - min_y)

            # Label should be benchmark name, not tag.  Really only need one 
            # legend, preferably off to the side.

            axis.set_title(tag)
            axis.scatter(data.rmsd, data.score,
                    color=colors.from_cycle(0), marker='o', edgecolor='none',
                    label=tag)

            axis.set_xlim(0, 10)
            axis.set_ylim(min_y, max_y)
            axis.axvline(1, color='gray', linestyle='--')


if __name__ == '__main__':
    import docopt
    arguments = docopt.docopt(__doc__)
    benchmarks = [int(x) for x in arguments['<benchmark_names_or_ids>']]

    data = extract_benchmark_data()

    plot_score_vs_rmsd_funnels(benchmarks, data)

    plt.show()

