#!/usr/bin/env python2

""" 
Generate a report comparing any number of benchmark runs.  The report will 
include the following figures for each benchmark being compared:

1. The distribution of the RMSD of the lowest scoring decoys.
2. The distribution of the fraction of sub-angstrom decoys sampled.
3. Score vs. RMSD plots for each structure.

Usage:
    analysis.py <benchmark_names_or_ids>...
"""

import re
import pandas as pd
import matplotlib.pyplot as plt

from helpers import settings; settings.load()
from helpers import database

if __name__ == '__main__':
    import docopt
    arguments = docopt.docopt(__doc__)
    benchmark_id = int(arguments['<benchmark_names_or_ids>'][0])

    with database.connect() as session:
        benchmarks = database.data_frame_from_table(session, database.Benchmarks)
        benchmark_protocols = database.data_frame_from_table(session, database.BenchmarkProtocols)
        protocols = database.data_frame_from_table(session, database.Protocols)
        batches = database.data_frame_from_table(session, database.Batches)
        structures = database.data_frame_from_table(session, database.Structures)
        rmsds = database.data_frame_from_table(session, database.ProteinRmsdNoSuperposition)
        scores = database.data_frame_from_table(session, database.TotalScores)

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
    scores_vs_rmsds = scores_vs_rmsds[
            scores_vs_rmsds.benchmark_id == benchmark_id]

    # Now I've got score vs rmsd for every structure...
    # make function for each plot
    # function loops through each benchmarks, caluclates metric, adds to plot

    tag = '1srp'
    pdb = 'structures/{0}.pdb'.format(tag)
    data = scores_vs_rmsds[scores_vs_rmsds.pdb == pdb]
    plt.title(pdb)
    plt.scatter(data.rmsd, data.score)

    plt.savefig('{0}.pdf'.format(tag))
    plt.show()
