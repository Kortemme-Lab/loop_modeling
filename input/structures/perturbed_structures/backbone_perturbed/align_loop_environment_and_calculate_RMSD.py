#!/usr/bin/env python2
'''Align two structures by the surrounding residues to a defined loop.
Then calculate the RMSD for the input loop and the grafted loop.

Usage:
    ./align_loop_environment_and_calculate_lDDT.py input_pdb ref_pdb loops_file
'''

import sys

import numpy as np

import pyrosetta
from pyrosetta import rosetta

def residue_bb_distance(pose, res1, res2):
    '''Return the heavy atom backbone distance between two residues'''
    bb_atoms = ['N', 'CA', 'C', 'O']

    return min(pose.residue(res1).xyz(a).distance(pose.residue(res2).xyz(a))
            for a in bb_atoms)

def get_loop_surrounding_residues(pose, loop_residues, cutoff_distance=10):
    '''Get residues around a loop.'''
    surrounding_residues = set()
    for i in range(1, pose.size() + 1):
        if i in loop_residues: continue

        for j in loop_residues:
            if residue_bb_distance(pose, i, j) < cutoff_distance:
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

def align_points_and_calculate_RMSD(xyzs1, xyzs2):
    '''Align the points and calculate the RMSD after alignment.'''
    def xyzV_to_np_array(xyz):
        return np.array([xyz.x, xyz.y, xyz.z])

    def RMSD(P1, P2):
        diff = [P1[i] - P2[i] for i in range(len(P1))]
        return np.sqrt(sum(np.dot(d, d) for d in diff) / len(P1))

    P1 = [xyzV_to_np_array(xyz) for xyz in xyzs1]
    P2 = [xyzV_to_np_array(xyz) for xyz in xyzs2]

    M, v = get_superimpose_transformation(P1, P2)

    P1_aligned = [np.dot(M, p) + v for p in P1]

    return RMSD(P1_aligned, P2)

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

def calc_bb_rmsd(pose1, pose2, residues):
    '''Calculate the backbone RMSD between two poses.'''
    bb_atoms = ['N', 'CA', 'C', 'O']

    dss = []
    for res in residues:
        for a in bb_atoms:
            dss.append(pose1.residue(res).xyz(a).distance_squared(pose2.residue(res).xyz(a)))

    return np.sqrt(sum(dss) / len(dss))

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
    surrounding_residues = [r for r in surrounding_residues_ref if r < pdb_pose.size()]

    align_poses(pdb_pose, ref_pose, surrounding_residues)

    print "loop rmsd after alignment = {0:.2f}".format(calc_bb_rmsd(pdb_pose, ref_pose, loop_residues))

    #pdb_pose.dump_pdb('test_aligned.pdb')###DEBUG
    #ref_pose.dump_pdb('test_ref.pdb')###DEBUG

    #pdb_loop_cas = [pdb_pose.residue(i).xyz('CA') for i in loop_residues]
    #ref_loop_cas = [ref_pose.residue(i).xyz('CA') for i in loop_residues]
    #pdb_env_cas = [pdb_pose.residue(i).xyz('CA') for i in surrounding_residues]
    #ref_env_cas = [ref_pose.residue(i).xyz('CA') for i in surrounding_residues]

    #print align_points_and_calculate_RMSD(pdb_env_cas, ref_env_cas)
    #print align_points_and_calculate_RMSD(pdb_loop_cas + pdb_env_cas, ref_loop_cas + ref_env_cas)
    #print align_points_and_calculate_RMSD(ref_loop_cas + pdb_env_cas, ref_loop_cas + ref_env_cas)
