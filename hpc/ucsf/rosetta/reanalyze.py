#!/usr/bin/env python2
'''Re-analyze the data of a data set. For example, 
re-calculate the RMSD of the loops using a different
aligning method. Usage:
    ./reanalyze.py benchmark_id reanalyze_method
'''

import os
import sys

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

def calc_rmsd_for_one_model(benchmark_id, model, reanalyze_method):
    '''Calculate the RMSD for one model defined as a tuple
    (pdb_id, model_id, loop_rmsd, total_energy, runtime)
    '''
    data_controller = DataController('disk')
    benchmark_define_dict = data_controller.get_benchmark_define_dict(benchmark_id)
    
    # Get the loops definition and the reference structure path
    
    for i, b_input in enumerate(benchmark_define_dict['input_pdbs']):
        if model[0] == os.path.basename(b_input.pdb_path)[:4]:
            ref_path = os.path.join(os.path.dirname(os.path.dirname(b_input.pdb_path)),
                    'reference', model[0] + '.pdb')
            protein_id = i
            loop_file = b_input.pdb_path[:-4] + '.loop'
            break

    # Get the model structure
    
    num_total_proteins = len(benchmark_define_dict['input_pdbs'])
    global_model_id = num_total_proteins * model[1] + protein_id
    model_path = os.path.join('data', benchmark_id, 'structures', '{0}_{1}_0001.pdb.gz'.format(global_model_id, model[0]))

    # Calculate RMSD
    
    #print ref_path, protein_id, model_path, loop_file###DEBUG
    print data_controller.calc_rmsd(loop_file, ref_path, model_path, rmsd_calculation_method=reanalyze_method)

if __name__ == '__main__':
    
    benchmark_id = sys.argv[1]
    reanalyze_method = sys.argv[2]

    
    results = load_existing_results(benchmark_id)

    calc_rmsd_for_one_model(benchmark_id, results[0], reanalyze_method)
