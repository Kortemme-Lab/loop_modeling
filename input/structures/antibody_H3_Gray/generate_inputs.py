#!/usr/bin/env python2

import os
import gzip
import subprocess
import re

import pyrosetta
from pyrosetta import rosetta

structures = [
        ('2ddq', ['H', 'L']),
        ('1dqq', ['A', 'B']),
        ('1z3g', ['H', 'L']),
        #('2ai0', ['H', 'L']),#Excluded by sphinx
        ('1tet', ['H', 'L']),
        ('1bql', ['H', 'L']),
        ('1cgs', ['H', 'L']),
        ('1mlb', ['A', 'B']),
        #('2cip', ['H', 'L']),#Not an antibody
        ('2bdn', ['H', 'L']),
        ('1fgn', ['H', 'L']),
        ('1jpt', ['H', 'L']),
        ('1a6t', ['A', 'B']),
        ('1kem', ['H', 'L']),
        ('1qbl', ['H', 'L']),
        ('1vfa', ['A', 'B']),
        ('1iqd', ['A', 'B']),
        ('1k4c', ['A', 'B']),
        ('1jhl', ['H', 'L']),
        ('2aep', ['H', 'L']),
        ('2fbj', ['H', 'L']),
        ('1igt', ['A', 'B']),
        ('2fd6', ['H', 'L']),
        ('2adf', ['H', 'L']),
        ('2jel', ['H', 'L']),
        #('1ynt', ['H', 'L']),#Rosetta cannot read Error with SSBond record.
        ('1dba', ['H', 'L']),
        ('2b2x', ['H', 'L']),
        ('1clz', ['H', 'L']),
        ('2cju', ['H', 'L']),
        ('1for', ['H', 'L']),
        ('1kb5', ['H', 'L']),
        ('2aju', ['H', 'L']),
        ('1ztx', ['H', 'L']),
        ('1mcp', ['H', 'L']),
        ('2fjg', ['H', 'L']),
        ('2h1p', ['H', 'L']),
        ('1fpt', ['H', 'L']),
        ('2adg', ['A', 'B']),
        ('1igm', ['H', 'L']),
        ('2fjh', ['H', 'L']),
        ('2g5b', ['H', 'L']),
        ('2h2p', ['C', 'D']),
        ('1fbi', ['H', 'L']),
        ('1bj1', ['H', 'L']),
        ('1wc7', ['H', 'L']),
        ('1zan', ['H', 'L']),
        ('2aj3', ['A', 'B']),
        ('2dqu', ['H', 'L']),
        ('1f58', ['H', 'L']),
        ('1hzh', ['H', 'L']),
        ('1g9m', ['H', 'L']),
        ('2b4c', ['H', 'L']),
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

def find_H3_loop(pose):
    '''Return the start and stop (inclusive) positions of
    the H3 loop in a pose.'''
    sequence = pose.sequence()

    query_seq = r'C.{2}(.{3,25})WG.G'
    matches = [m for m in re.finditer(query_seq, sequence)]
    
    print matches
    assert(len(matches) == 1)

    m = matches[0]
    
    loop = (m.start() + 4, m.start() + 3 + len(m.group(1)))
    return loop

def prepare_inputs(pdb_id, chains_to_keep):
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
    pose.pdb_info().obsolete(True)

    # Get the loop

    print pdb_id
    loop = find_H3_loop(pose)
    
    with open(os.path.join('inputs', pdb_id + '.loop'), 'w') as f:
        f.write('LOOP {0} {1} {1} 0 1'.format(loop[0], loop[1]))

    # Generate the fasta file, this is important because the chains to be modeled are not labeled as chain A

    with open(os.path.join('inputs', pdb_id + '.fasta'), 'w') as f:
        f.write('>{0}|A\n'.format(pdb_id))
        f.write(pose.sequence())

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
        prepare_inputs(structure[0], structure[1])
