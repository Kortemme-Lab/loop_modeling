#!/usr/bin/env python2

import os
import gzip
import subprocess
import re

import pyrosetta
from pyrosetta import rosetta

structures = [
        ('4ma3', ['H', 'L'], 'ARGLYNDYTV'),
	('4kuz', ['H', 'L'], 'ARNDGYRGYAMDY'),	
	#('4kq3', ['H', 'L'], 'ARHWGQGTLVTVSSH'),	
	('4kq3', ['H', 'L'], 'ARDSEYYFDH'),	
	('4kq4', ['H', 'L'], 'AREVRRSMDY'),
	('4m6m', ['H', 'L'], 'ARWYYKPFDV'),	
	('4m6o', ['H', 'L'], 'ARVYSSGWHVSDYFDY'),
	('4mau', ['H', 'L'], 'ASYDGYSFDY'),
	('4m7k', ['H', 'L'], 'ARSGYYGNSGFAY'),	
	('4kmt', ['H', 'L'], 'ARYDGIYGELDF'),
	('4m61', ['A', 'B'], 'ARGRLRRGGYFDY'),
	('4m43', ['H', 'L'], 'ARRNYDGSWFAY')
        ]

def download_pdb(pdb_id):
    '''Download a pdb file'''
    raw_file = 'pdb{0}.ent.gz'.format(pdb_id)
    url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/{0}/{1}'.format(pdb_id[1:3], raw_file)

    subprocess.call(['wget', url])

    with gzip.open(raw_file, 'r') as fin:
        with open(pdb_id + '.pdb', 'w') as fout:
            fout.write(fin.read())

    os.remove(raw_file)

def find_sequence_in_pose(pose, query_seq):
    '''Return the start and stop (inclusive) positions of
    a query sequence in a pose.'''
    sequence = pose.sequence()

    positions = [m.start() for m in re.finditer(query_seq, sequence)]

    print positions
    assert(len(positions) == 1)

    return positions[0] + 1, positions[0] + len(query_seq)


def prepare_inputs(pdb_id, chains_to_keep, h3_sequence):
    '''Prepare for the input files.'''
    download_pdb(pdb_id)

    # Extract the chains to keep

    raw_pose = rosetta.core.pose.Pose()
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(raw_pose, pdb_id + '.pdb')
   
    seqposes_to_keep = rosetta.utility.vector1_unsigned_long()
    for seqpos in range(1, raw_pose.size()):
        if raw_pose.pdb_info().chain(seqpos) in chains_to_keep:
            seqposes_to_keep.append(seqpos)
    
    rosetta.core.pose.pdbslice(pose, raw_pose, seqposes_to_keep)
    pose.pdb_info().obsolete()

    # Get the loop

    print pdb_id
    loop = find_sequence_in_pose(pose, h3_sequence)
    
    with open(os.path.join('inputs', pdb_id + '.loop'), 'w') as f:
        f.write('LOOP {0} {1} {1} 0 1'.format(loop[0], loop[1]))

    pose.dump_file(os.path.join('inputs', pdb_id + '.pdb'))
    pose.dump_file(os.path.join('reference', pdb_id + '.pdb'))

    os.remove(pdb_id + '.pdb')

if __name__ == '__main__':
    pyrosetta.init(options='-ignore_unrecognized_res')

    if not os.path.exists('inputs'):
        os.mkdir('inputs')
    if not os.path.exists('reference'):
        os.mkdir('reference')

    for structure in structures:
        prepare_inputs(structure[0], structure[1], structure[2])
