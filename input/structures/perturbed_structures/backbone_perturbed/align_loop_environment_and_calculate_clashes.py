#!/usr/bin/env python2
'''Align two structures by the surrounding residues to a defined loop.
Then calculate clashes for the input loop and the grafted loop.

Usage:
    ./align_loop_environment_and_calculate_lDDT.py input_pdb ref_pdb loops_file
'''

import sys

import numpy as np

import pyrosetta
from pyrosetta import rosetta

def residue_bb_distance(res1, res2):
    '''Return the heavy atom backbone distance between two residues'''
    bb_atoms = ['N', 'CA', 'C', 'O']

    return min(res1.xyz(a).distance(res2.xyz(a))
            for a in bb_atoms)

def get_loop_surrounding_residues(pose, loop_residues, cutoff_distance=10):
    '''Get residues around a loop.'''
    surrounding_residues = set()
    for i in range(1, pose.size() + 1):
        if i in loop_residues: continue

        for j in loop_residues:
            if residue_bb_distance(pose.residue(i), pose.residue(j)) < cutoff_distance:
                surrounding_residues.add(i)

    return list(surrounding_residues)

def get_superimpose_transformation(P1, P2):
    '''Get the superimpose transformation that transfoms a list of
    points P1 to another list of points P2.'''
    if len(P1) != len(P2):
        raise Exception("Sets to be superimposed must have same number of points.")

    com1 = np.mean(P1, axis=0)
    com2 = np.mean(P2, axis=0)

    R = np.dot(np.transpose(np.array(P1) - com1), np.array(P2) - com2)
    V, S, W = np.linalg.svd(R)

    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]

    M = np.transpose(np.array(np.dot(V, W)))

    return M, com2 - np.dot(M, com1)

def align_poses(source_pose, target_pose, residues_to_align):
    '''Align the source pose to the target pose.
    Only backbone heavey atoms are aligned.
    '''
    def xyzV_to_np_array(xyz):
        return np.array([xyz.x, xyz.y, xyz.z])

    def np_array_to_xyzV(a):
        return rosetta.numeric.xyzVector_double_t(a[0], a[1], a[2])
    
    def xyzM_to_np_array(M):
        return np.array([[M.xx, M.xy, M.xz],
                         [M.yx, M.yy, M.yz],
                         [M.zx, M.zy, M.zz]])
    
    def np_array_to_xyzM(a):
        return rosetta.numeric.xyzMatrix_double_t.rows(
                a[0][0], a[0][1], a[0][2],
                a[1][0], a[1][1], a[1][2],
                a[2][0], a[2][1], a[2][2])

    bb_atoms = ['N', 'CA', 'C', 'O']

    points_source = [xyzV_to_np_array(source_pose.residue(r).xyz(a))
            for r in residues_to_align for a in bb_atoms]
    points_target = [xyzV_to_np_array(target_pose.residue(r).xyz(a))
            for r in residues_to_align for a in bb_atoms]

    R, v = get_superimpose_transformation(points_source, points_target)

    source_pose.apply_transform_Rx_plus_v(np_array_to_xyzM(R), np_array_to_xyzV(v))

def find_residue_bb_clashes(residues, cutoff_distance=2.5):
    '''Find backbone clashes within a list of residues.'''
    N = len(residues)

    clahes_detected = []

    for i in range(N):
        for j in range(i + 1, N):
            if -1 <= residues[i].seqpos() - residues[j].seqpos() <= 1: continue

            d = residue_bb_distance(residues[i], residues[j]) 
            if d < cutoff_distance:
                clahes_detected.append((residues[i].seqpos(), residues[j].seqpos(), d))

    return clahes_detected

def residue_heavy_atom_clashes(residue1, residue2, cutoff_distance=2.5):
    '''Return true if two residues have heavy atoms that clash.'''
    for i in range(1, residue1.nheavyatoms() + 1):
        p1 = residue1.xyz(i)
        for j in range(1, residue2.nheavyatoms() + 1):
            p2 = residue2.xyz(j)

            if p1.distance(p2) < cutoff_distance:
                return True

    return False

def replace_intra_residue_torsions(pose, seqpos, ref_residue):
    '''Replace the intra residue torsions at seqpos
    by the torsions from the reference residue.
    '''
    ref_chis = ref_residue.chi()
    
    for i in range(1, pose.residue(seqpos).nchi() + 1):
        pose.set_chi(i, seqpos, ref_chis[i])

def residues_clash_unavoidable_clash(pose1, res1, pose2, res2):
    '''Return True if the clash between two residues is unavoidable.'''
    rotamer_set1 = rosetta.core.pack.rotamer_set.bb_independent_rotamers( pose1.residue(res1).type(), True )
    rotamer_set2 = rosetta.core.pack.rotamer_set.bb_independent_rotamers( pose2.residue(res2).type(), True )

    for i in range(len(rotamer_set1) + 1):
        # When i == 0, use the current rotamer

        if i > 0:
            replace_intra_residue_torsions(pose1, res1, rotamer_set1[i])

        for j in range(len(rotamer_set2) + 1):
            # When j == 0, use the current rotamer
            
            if j > 0:
                replace_intra_residue_torsions(pose2, res2, rotamer_set2[j])

            if not residue_heavy_atom_clashes(pose1.residue(res1), pose2.residue(res2)):
                return False

    return True

def find_unavoidable_heavy_atoms_clahses(pose1, residues1, pose2, residues2):
    '''Find unavoidable clashes between two list for residues.'''
    clashing_residues = []
    
    for res1 in residues1:
        for res2 in residues2:
            if -1 <= res1 - res2 <= 1: continue

            if residues_clash_unavoidable_clash(pose1, res1, pose2, res2):
                clashing_residues.append((res1, res2))

    return clashing_residues


if __name__ == '__main__':

    pyrosetta.init()

    pdb_file = sys.argv[1]
    ref_file = sys.argv[2]
    loops_file = sys.argv[3] 

    pdb_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pdb_pose, pdb_file)
    ref_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(ref_pose, ref_file)

    loops = rosetta.protocols.loops.Loops(loops_file)

    loop_residues = []
    for i in range(loops.size()):
        loop_residues += list(range(loops[i + 1].start(), loops[i + 1].stop() + 1))

    surrounding_residues_ref = get_loop_surrounding_residues(ref_pose, loop_residues)
    surrounding_residues = [r for r in surrounding_residues_ref if r <= pdb_pose.size()]

    align_poses(pdb_pose, ref_pose, surrounding_residues)

    pdb_pose.dump_pdb('test_aligned.pdb')###DEBUG
    ref_pose.dump_pdb('test_ref.pdb')###DEBUG

    #pdb_loop_res = [pdb_pose.residue(i) for i in loop_residues]
    #ref_loop_res = [ref_pose.residue(i) for i in loop_residues]
    #pdb_env_res = [pdb_pose.residue(i) for i in surrounding_residues]
    #ref_env_res = [ref_pose.residue(i) for i in surrounding_residues]

    #print 'BB clashes for the input loop', find_residue_bb_clashes(pdb_loop_res + pdb_env_res)
    #print 'BB clashes for the grafted loop', find_residue_bb_clashes(ref_loop_res + pdb_env_res)

    input_input_clashes = find_unavoidable_heavy_atoms_clahses(pdb_pose, loop_residues, pdb_pose, surrounding_residues)
    native_native_clashes = find_unavoidable_heavy_atoms_clahses(ref_pose, loop_residues, ref_pose, surrounding_residues)
    grafted_input_clashes = find_unavoidable_heavy_atoms_clahses(ref_pose, loop_residues, pdb_pose, surrounding_residues)
    
    print 'Unavoidable clashes between the input loop and the input environment are', input_input_clashes 
    print 'Unavoidable clashes between the native loop and the native environment are', native_native_clashes
    print 'Unavoidable clashes between the grafted loop and the input environment are', grafted_input_clashes
