import sys
import os
import math
import glob
import json
import time
import pprint

import pandas
import numpy

from tools import colortext
from tools.fs.fsio import read_file, get_file_lines, write_file
from tools.bio.pdb import PDB
from tools.bio.basics import backbone_atoms
from tools.pymath.cartesian.rmsd import compute_rmsd_by_matrix


def compute_rmsds(results_folder, expectn, top_x):

    assert(top_x <= expectn)

    # Set up reference structures
    structures_folder = '../input/structures'
    rcsb_references = os.path.join(structures_folder, 'rcsb', 'reference')
    rosetta_references = os.path.join(structures_folder, 'rosetta', 'reference')

    best_scoring_structures = {}
    median_scoring_structures = {}
    worst_scoring_structures = {}
    percent_subangstrom = {}
    top_x_structures = {}

    csv_file = ['\t'.join(['PDB ID', 'Models', '%<1.0A', 'Best score', 'Top{0} score'.format(top_x), 'Median score', 'Worst score', 'Closest score', 'Top1 RMSD', 'Top{0} RMSD'.format(top_x), 'Closest RMSD'])]
    pdb_ids = [os.path.splitext(os.path.split(s.strip())[1])[0] for s in get_file_lines('../input/full.pdbs') if s.strip()]
    for pdb_id in pdb_ids:

        colortext.message('Computing RMSDs for {0}.'.format(pdb_id))
        rcsb_reference_pdb = os.path.join(rcsb_references, pdb_id + '.pdb')
        assert(os.path.exists(rcsb_reference_pdb))
        rosetta_reference_pdb = os.path.join(rcsb_references, pdb_id + '.pdb')
        assert(os.path.exists(rosetta_reference_pdb))
        assert(len(pdb_id) == 4)
        loops_file = os.path.join('../input/structures/rosetta/pruned/{0}.loop.json'.format(pdb_id))
        loop_sets = json.loads(read_file(loops_file))
        assert(len(loop_sets['LoopSet']) == 1)

        # Parsing the score files and ordering the structures by score
        total_scores = {} # a mapping from scores to run_ids / structure_ids
        pdb_loop_residue_matrices = {}
        run_ids = []
        for sc_file in glob.glob(os.path.join(results_folder, '{0}*.sc'.format(pdb_id))):
            sc_filename = os.path.split(sc_file)[1]
            assert(sc_filename.startswith('{0}_score'.format(pdb_id)))
            run_id = int(sc_filename[10:-3])
            sc_lines = [l.strip() for l in get_file_lines(sc_file) if l.strip()]
            assert(sc_lines[0] == 'SEQUENCE:')
            assert(sc_lines[1].split()[:2] == ['SCORE:', 'total_score'])
            assert(sc_lines[2].split()[0] == 'SCORE:')
            total_score = float(sc_lines[2].split()[1])
            total_scores[total_score] = total_scores.get(total_score, [])
            total_scores[total_score].append(run_id)
            run_ids.append(run_id)

            # Extract the PDB coordinates into a pandas dataframe (HDF5 format)
            associated_pdb_file = os.path.join(results_folder, '{0}_{0}{1}_0001.pdb'.format(pdb_id, run_id))
            assert(os.path.exists(associated_pdb_file))
            hdf5_file = os.path.splitext(associated_pdb_file)[0] + '.hdf5'
            if os.path.exists(hdf5_file):
                store = pandas.HDFStore(hdf5_file)
                pdb_loop_residue_matrix = store['dataframe']
                store.close()
            else:
                pdb_loop_residue_matrix = PDB.extract_xyz_matrix_from_loop_json(PDB.from_filepath(associated_pdb_file).structure_lines, loop_sets, atoms_of_interest = backbone_atoms, expected_num_residues = 12, expected_num_residue_atoms = 4)
                store = pandas.HDFStore(hdf5_file)
                store['dataframe'] = pdb_loop_residue_matrix
                store.close()
            pdb_loop_residue_matrices[run_id] = pdb_loop_residue_matrix

        num_structures = len(run_ids)
        if num_structures < expectn:
            print('Error: Expected {0} structures but only found {1}.'.format(expectn, num_structures))
            sys.exit(1)

        # Truncate the structures to the top expectn-scoring files
        structure_ids_by_score = []
        structure_scores_by_id = {}
        for score, structure_ids in sorted(total_scores.iteritems()):
            for structure_id in sorted(structure_ids): # sorting here ensures determinism
                structure_ids_by_score.append((score, structure_id))
                structure_scores_by_id[structure_id] = score
        structure_ids_by_score = structure_ids_by_score[:expectn]
        assert(len(structure_ids_by_score) == expectn)

        # Determine the top-X-scoring structures and the median-scoring structure
        top_x_structures[pdb_id] = structure_ids_by_score[:top_x]
        median_scoring_structures[pdb_id] = structure_ids_by_score[int(expectn / 2)]

        # Determine the lowest-/best-scoring structure
        best_scoring_structures[pdb_id] = structure_ids_by_score[0]
        best_score = best_scoring_structures[pdb_id][0]
        print(top_x_structures[pdb_id])
        print(best_scoring_structures[pdb_id])
        assert(top_x_structures[pdb_id][0] == best_scoring_structures[pdb_id])

        # Determine the worst-scoring structure
        worst_scoring_structures[pdb_id] = structure_ids_by_score[-1]
        worst_score = worst_scoring_structures[pdb_id][0]

        print('best', best_scoring_structures[pdb_id])
        print('median', median_scoring_structures[pdb_id])
        print('worst', worst_scoring_structures[pdb_id])

        # Read the coordinates from the reference PDB file
        rcsb_reference_matrix = PDB.extract_xyz_matrix_from_loop_json(PDB.from_filepath(rcsb_reference_pdb).structure_lines, loop_sets, atoms_of_interest = backbone_atoms, expected_num_residues = 12, expected_num_residue_atoms = 4)
        rosetta_reference_matrix = PDB.extract_xyz_matrix_from_loop_json(PDB.from_filepath(rosetta_reference_pdb).structure_lines, loop_sets, atoms_of_interest = backbone_atoms, expected_num_residues = 12, expected_num_residue_atoms = 4)

        # Compute RMSDs for all expectn structures using the pandas dataframes
        c = 0
        num_subangstrom = 0
        structure_rmsds = {}
        for score, structure_ids in sorted(total_scores.iteritems()):
            for structure_id in structure_ids:
                predicted_matrix = pdb_loop_residue_matrices[structure_id]
                structure_rmsds[structure_id] = compute_rmsd_by_matrix(rcsb_reference_matrix, predicted_matrix)
                if structure_rmsds[structure_id] < 1.0:
                    num_subangstrom += 1
                assert(structure_rmsds[structure_id] == compute_rmsd_by_matrix(rosetta_reference_matrix, predicted_matrix))
        assert(len(structure_rmsds) == expectn)
        percent_subangstrom[pdb_id] = (float(num_subangstrom) / float(expectn)) * 100.0

        closest_score = worst_score + 1.0
        closest_rmsd = min(structure_rmsds.values())
        for structure_id, structure_rmsd in structure_rmsds.iteritems():
            if structure_rmsd == closest_rmsd:
                closest_score = min(closest_score, structure_scores_by_id[structure_id])

        #best_scoring_pdb_coordinate_matrix = pdb_loop_residue_matrices[best_scoring_structures[pdb_id][1]]
        #print(best_scoring_pdb_coordinate_matrix)
        #sys.exit(0)

        #best_scoring_pdb_file = os.path.join(results_folder, '{0}_{0}{1}_0001.pdb'.format(pdb_id, best_scoring_structures[pdb_id][1]))
        #print(best_scoring_pdb_file)

        #assert(os.path.exists(best_scoring_pdb_file))

        # Read the start and stop residue IDs from file
        #start_pdb_residue_id = PDB.ChainResidueID2String(loop_sets['LoopSet'][0]['start']['chainID'], str(loop_sets['LoopSet'][0]['start']['resSeq']) + loop_sets['LoopSet'][0]['start']['iCode'])
        #stop_pdb_residue_id = PDB.ChainResidueID2String(loop_sets['LoopSet'][0]['stop']['chainID'], str(loop_sets['LoopSet'][0]['stop']['resSeq']) + loop_sets['LoopSet'][0]['stop']['iCode'])



        # Compute RMSDs for all Top expectn structures using Python
        #python_rmsd_time = time.time()
        #reference_matrix = PDB.extract_xyz_matrix_from_pdb_residue_range(PDB.from_filepath(rcsb_reference_pdb).structure_lines, start_pdb_residue_id, stop_pdb_residue_id, atoms_of_interest = backbone_atoms, expected_num_residues = 12, expected_num_residue_atoms = 4)
        #c = 0
        #for sc, structure_ids in sorted(total_scores.iteritems()):
        #    for structure_id in structure_ids:
        #        c += 1
        #        if c <= expectn:
        #            associated_pdb_file = os.path.join(results_folder, '{0}_{0}{1}_0001.pdb'.format(pdb_id, structure_id))
        #            rcsb_rmsd = calculate_rmsd_from_pdb_files(associated_pdb_file, rcsb_reference_pdb, start_pdb_residue_id, stop_pdb_residue_id, expected_num_residues = 12, expect_all_atoms_per_residue = True)
        #            #print(c, rcsb_rmsd)
        #python_rmsd_time = time.time() - python_rmsd_time
        #print('RMSD-by-python computation time for {0} structures: {1}s.'.format(expectn, python_rmsd_time))
        #sys.exit(0)

        #predicted_matrix = PDB.extract_xyz_matrix_from_pdb_residue_range(PDB.from_filepath(best_scoring_pdb_file).structure_lines, start_pdb_residue_id, stop_pdb_residue_id, atoms_of_interest = backbone_atoms, expected_num_residues = 12, expected_num_residue_atoms = 4)

        #reference_matrix_2 = PDB.extract_xyz_matrix_from_loop_json(PDB.from_filepath(rcsb_reference_pdb).structure_lines, loop_sets, atoms_of_interest = backbone_atoms, expected_num_residues = 12, expected_num_residue_atoms = 4)
        #predicted_matrix_2 = PDB.extract_xyz_matrix_from_loop_json(PDB.from_filepath(best_scoring_pdb_file).structure_lines, loop_sets, atoms_of_interest = backbone_atoms, expected_num_residues = 12, expected_num_residue_atoms = 4)

        top_1_rmsd = structure_rmsds[best_scoring_structures[pdb_id][1]]
        top_x_rmsd = min([structure_rmsds[top_x_structures[pdb_id][j][1]] for j in range(top_x)])
        top_x_score = worst_score + 1.0
        for j in range(top_x):
            top_x_score = min(top_x_score, top_x_structures[pdb_id][j][0])
        assert(top_x_score <= worst_score)
        assert(top_x_rmsd <= top_1_rmsd)

        #rcsb_rmsd = calculate_rmsd_from_pdb_files(best_scoring_pdb_file, rcsb_reference_pdb, start_pdb_residue_id, stop_pdb_residue_id, expected_num_residues = 12, expect_all_atoms_per_residue = True)
        #rosetta_rmsd = calculate_rmsd_from_pdb_files(best_scoring_pdb_file, rosetta_reference_pdb, start_pdb_residue_id, stop_pdb_residue_id, expected_num_residues = 12, expect_all_atoms_per_residue = True)
        #rcsb_rmsd_2 = compute_rmsd_by_matrix(reference_matrix, predicted_matrix)
        #rcsb_rmsd_3 = compute_rmsd_by_matrix(reference_matrix_2, predicted_matrix_2)
        #print(rcsb_rmsd, rosetta_rmsd, rcsb_rmsd_2)
        #assert(rcsb_rmsd == rosetta_rmsd and abs(rcsb_rmsd - rcsb_rmsd_2) < 0.0001 and abs(rcsb_rmsd - rcsb_rmsd_3) < 0.0001)
        print('Top 1 RMSD (predicted vs Rosetta/RCSB reference structure): {0}'.format(top_1_rmsd))
        print('Top {0} RMSD (predicted vs Rosetta/RCSB reference structure): {1}'.format(top_x, top_x_rmsd))

        #top_5_rmsds = []
        #for structure_id in best_five_scoring_structures:
        #    pdb_file = os.path.join(results_folder, '{0}_{0}{1}_0001.pdb'.format(pdb_id, structure_id))
        #    top_5_rmsds.append(calculate_rmsd_from_pdb_files(pdb_file, rcsb_reference_pdb, start_pdb_residue_id, stop_pdb_residue_id))
        #top_5_rmsd = min(top_5_rmsds)
        #print('Top5 RMSD vs Rosetta/RCSB: {0}'.format(top_5_rmsd))
        csv_file.append('\t'.join(map(str, [pdb_id, expectn, percent_subangstrom[pdb_id], best_score, top_x_score, median_scoring_structures[pdb_id][0], worst_score, closest_score, top_1_rmsd, top_x_rmsd, closest_rmsd])))

    write_file('analysis.csv', '\n'.join(csv_file))


def get_atoms(pdb_filepath, start_pdb_residue_id, stop_pdb_residue_id, atoms_of_interest = backbone_atoms, expected_num_residues = None, expected_num_residue_atoms = None):
    atoms = {}
    found_start = False
    found_end = False
    for l in PDB.from_filepath(pdb_filepath).structure_lines:
        res_id = None
        atom_type = None
        if l.startswith('ATOM  '):
            res_id = l[21:27]
            atom_type = l[12:16].strip()
            if res_id == start_pdb_residue_id:
                found_start = True
            if res_id == stop_pdb_residue_id:
                assert(found_start)
                found_end = True

        if found_end and res_id != stop_pdb_residue_id:
            break

        if found_start and l.startswith('ATOM  ') and (not(atoms_of_interest) or (atom_type in atoms_of_interest)):
            assert(res_id and atom_type and not(atoms.get(res_id, {}).get(atom_type)))
            atoms[res_id] = atoms.get(res_id, {})
            atoms[res_id][atom_type] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))

    if expected_num_residues != None:
        assert(len(atoms) == expected_num_residues)
    if expected_num_residue_atoms != None:
        for res_id, atom_details in atoms.iteritems():
            assert(len(atom_details) == expected_num_residue_atoms)
    return atoms


def calculate_rmsd_from_pdb_files(pdb1_filepath, pdb2_filepath, start_pdb_residue_id, stop_pdb_residue_id, atoms_of_interest = backbone_atoms, expected_num_residues = None, expect_all_atoms_per_residue = False):

    expected_num_residue_atoms = None
    if expect_all_atoms_per_residue:
        expected_num_residue_atoms = len(atoms_of_interest)

    atoms1 = get_atoms(pdb1_filepath, start_pdb_residue_id, stop_pdb_residue_id, atoms_of_interest = atoms_of_interest, expected_num_residues = expected_num_residues, expected_num_residue_atoms = expected_num_residue_atoms)
    atoms2 = get_atoms(pdb2_filepath, start_pdb_residue_id, stop_pdb_residue_id, atoms_of_interest = atoms_of_interest, expected_num_residues = expected_num_residues, expected_num_residue_atoms = expected_num_residue_atoms)

    for k, v in atoms1.iteritems():
        if expect_all_atoms_per_residue:
            assert(len(v) == len(atoms_of_interest))
    for k, v in atoms2.iteritems():
        if expect_all_atoms_per_residue:
            assert(len(v) == len(atoms_of_interest))
    assert(atoms1.keys() == atoms2.keys())

    c = 0.0
    total = 0.0
    for res_id, atoms in sorted(atoms1.iteritems()):
        for atom, xyz1 in sorted(atoms.iteritems()):
            xyz2 = atoms2[res_id].get(atom)
            if expect_all_atoms_per_residue:
                assert(xyz2)
            if xyz2:
                c += 1.0
                for w in range(3):
                    total += ((xyz2[w] - xyz1[w]) * (xyz2[w] - xyz1[w]))
    total = float(total) / float(c)
    total = math.sqrt(total)
    return (total)

if __name__ == '__main__':
    # todo: replace with docopt
    expectn = 500
    if len(sys.argv) == 3:
        try:
            expectn = int(sys.argv[2])
        except:
            print('Error: The second argument must be the expected number of structures for case.')
            sys.exit(1)
        sys.argv = sys.argv[:2]
    if len(sys.argv) == 2:
        structure_path = sys.argv[1]
        if not os.path.exists(structure_path):
            print('Error: The path "{0}" does not exist.'.format(structure_path))
            sys.exit(1)
        compute_rmsds(structure_path, expectn, 5)
    else:
        print('Error: Please specify a path containing the .sc and .pdb files e.g. "{0} /some/path".'.format(sys.argv[0]))
        sys.exit(1)
