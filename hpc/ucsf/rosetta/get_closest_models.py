#!/usr/bin/env python2
'''Get the closest models of a data set. Save
the closest models into a new directory in that dataset
    ./reanalyze.py benchmark_id
'''

import os
import sys; sys.path.append(os.getcwd())
import shutil

from libraries.dataController import DataController 

def load_existing_results(benchmark_id):
    '''Load the result file into a list of tuples
    (pdb_id, model_id, loop_rmsd, total_energy, runtime)
    '''
    # Find the result file

    for f in os.listdir(os.path.join('data', benchmark_id)):
        if f.endswith('.results'):
            result_file = os.path.join('data', benchmark_id, f)
            break

    # Load the result file

    results = []
    with open(result_file, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i > 0:
                sl = line.split()
                results.append((sl[0], int(sl[1]), float(sl[2]), float(sl[3]), int(sl[4])))

    return results 

if __name__ == '__main__':
    
    benchmark_id = sys.argv[1]

    # Load the existing results

    results = load_existing_results(benchmark_id)

    data_controller = DataController('disk')
    benchmark_define_dict = data_controller.get_benchmark_define_dict(benchmark_id)

    pdb_to_ids = {}
    pdb_to_lowest_rmsd = {}
    pdb_to_lowest_rmsd_model_id = {}

    for i, b_input in enumerate(benchmark_define_dict['input_pdbs']):
        pdb = os.path.basename(b_input.pdb_path)[:4]
        pdb_to_ids[pdb] = i
        pdb_to_lowest_rmsd[pdb] = float('inf')
        pdb_to_lowest_rmsd_model_id[pdb] = 0

    # Find the lowest scoring models

    for model in results:
        pdb, model_id, rmsd = model[0], model[1], model[2]
        
        if rmsd < pdb_to_lowest_rmsd[pdb]:
            pdb_to_lowest_rmsd[pdb] = rmsd 
            pdb_to_lowest_rmsd_model_id[pdb] = model_id

    # Save the lowest energy model to a new directory

    lowest_rmsd_model_dir = os.path.join('data', benchmark_id, 'lowest_rmsd_models')
    if not os.path.exists(lowest_rmsd_model_dir):
        os.mkdir(lowest_rmsd_model_dir)

    num_total_proteins = len(pdb_to_ids.keys())

    for pdb in pdb_to_ids.keys():
        global_model_id = num_total_proteins * pdb_to_lowest_rmsd_model_id[pdb] + pdb_to_ids[pdb]
       
        shutil.copy(os.path.join('data', benchmark_id, 'structures', '{0}_{1}_0001.pdb.gz'.format(global_model_id, pdb)),
                lowest_rmsd_model_dir)
