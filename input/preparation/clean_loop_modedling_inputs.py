#!/usr/bin/env python2

import os
import numpy as np

import pyrosetta
from pyrosetta import rosetta

def xyzV_to_np_array(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

def np_array_to_xyzV(a):
    return rosetta.numeric.xyzVector_double_t(a[0], a[1], a[2])

def np_array_to_xyzM(a):
    return rosetta.numeric.xyzMatrix_double_t.rows(
            a[0][0], a[0][1], a[0][2],
            a[1][0], a[1][1], a[1][2],
            a[2][0], a[2][1], a[2][2])

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

def get_backbone_points(pose, residues):
    '''Get backbone points for residues in a pose.'''
    points = []

    for res in residues:
        for atom in ['N', 'CA', 'C']:
            points.append(xyzV_to_np_array(pose.residue(res).xyz(atom)))

    return points

def superimpose_poses_by_residues(pose_source, residues_source, pose_target, residues_target):
    '''Superimpose residues in a source pose into residues in a target pose.
    Only backbone atoms are used for the superimposition.
    '''
    assert(len(residues_source) == len(residues_target))

    # Get the points to be superimposed

    points_source = get_backbone_points(pose_source, residues_source)
    points_target = get_backbone_points(pose_target, residues_target)

    # Get the rigid body transformation


    M, t = get_superimpose_transformation(points_source, points_target)

    # Transform the source pose

    pose_source.apply_transform_Rx_plus_v(np_array_to_xyzM(M), 
            np_array_to_xyzV(t))

def superimpose_reference_to_input(input_pdb_file, ref_pdb_file, loop_start_id, loop_stop_id):
    '''Superimpose the reference structure to the input structure by
    aligning the non-loop residues. Return the pose of the reference structure.
    '''
    pose_input = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose_input, input_pdb_file)
    pose_ref = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose_ref, ref_pdb_file)
   
    input_chain = pose_input.pdb_info().chain(1)
    input_start_seqpos = pose_input.pdb_info().pdb2pose(input_chain, loop_start_id[1])
    input_stop_seqpos = pose_input.pdb_info().pdb2pose(input_chain, loop_stop_id[1])
    ref_chain = pose_ref.pdb_info().chain(1)
    ref_start_seqpos = pose_ref.pdb_info().pdb2pose(ref_chain, loop_start_id[1])
    ref_stop_seqpos = pose_ref.pdb_info().pdb2pose(ref_chain, loop_stop_id[1])

    residues_input = []
    residues_ref = []

    for i in range(1, min(pose_input.size(), pose_ref.size()) + 1):
        if i < input_start_seqpos or i > input_stop_seqpos:
            residues_input.append(i)
        #if i < ref_start_seqpos or i > ref_stop_seqpos:
            residues_ref.append(i)

    superimpose_poses_by_residues(pose_ref, residues_ref, pose_input, residues_input)

    return pose_ref

def get_rosetta_loop_from_pdb_indices(pdb_file, start_id, stop_id, output_pdb_file, output_loop_file):
    '''Get the rosetta loop from the pdb indices. The start
    and stop indices are represented as pairs of (chain, atom_id).
    Write the structure and loop definition into new files.
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    chain = pose.pdb_info().chain(1)
    start_seqpos = pose.pdb_info().pdb2pose(chain, start_id[1])
    stop_seqpos = pose.pdb_info().pdb2pose(chain, stop_id[1])

    with open(output_loop_file, 'w') as f:
        f.write('LOOP {0} {1} {1} 0 1'.format(start_seqpos, stop_seqpos))

    pose.dump_pdb(output_pdb_file)


if __name__ == '__main__':
    pyrosetta.init(options='-ignore_unrecognized_res true')

    if not os.path.exists('inputs'):
        os.mkdir('inputs')
    if not os.path.exists('reference'):
        os.mkdir('reference')

    # Load the loops

    pdb_loops = []

    with open('trglist_8res', 'r') as f:
        for line in f.readlines():
            sl = line.split()
            start = sl[1].split(':')
            stop = sl[2].split(':')

            if start[0] == '_': start[0] = ' '
            if stop[0] == '_': stop[0] = ' '

            pdb_loops.append((sl[0], (start[0], int(start[1])), (stop[0], int(stop[1]))))

    # Move the native structures into the reference folder

    for pl in pdb_loops:
        print pl
        pose = superimpose_reference_to_input(pl[0] + '.bbpert.pdb', pl[0] + '.pdb', pl[1], pl[2])
        pose.dump_pdb(os.path.join('reference', pl[0] + '.pdb'))

    # Generate loop files and copy them into the inputs folder
    
    for pl in pdb_loops:
        get_rosetta_loop_from_pdb_indices(pl[0] + '.bbpert.pdb', pl[1], pl[2], os.path.join('inputs', pl[0] + '.pdb'), os.path.join('inputs', pl[0] + '.loop'))


