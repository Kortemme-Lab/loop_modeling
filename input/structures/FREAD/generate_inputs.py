#!/usr/bin/env python2

import os
import gzip
import subprocess
import re

import pyrosetta
from pyrosetta import rosetta

structures_12res = [
        ('1w5r', 'A', 12, 96), 
        ('1ae9', 'A', 12, 297),
        ('2i47', 'B', 12, 264),
        ('1jdl', 'A', 12, 93),
        ('2dxc', 'J', 12, 46), 
        ('1w96', 'B', 12, 493),
        ('1xcr', 'B', 12, 219),
        ('2iht', 'B', 12, 350),
        ('2o04', 'A', 12, 17),
        ('1eqc', 'A', 12, 20), 
        ('1yq2', 'E', 12, 335),
        ('1oew', 'A', 12, 171),
        ('1e5k', 'A', 12, 13),
        ('2cvd', 'D', 12, 136),
        ('1v7w', 'A', 12, 84), 
        ('1gyo', 'B', 12, 11), 
        ('2g8j', 'A', 12, 53), 
        ('1n1t', 'A', 12, 592),
        ('3bix', 'B', 12, 570),
        ('1qpc', 'A', 12, 455),
        ('2ywi', 'B', 12, 172),
        ('1rcq', 'A', 12, 157),
        ('2bo4', 'B', 12, 211),
        ('1dj0', 'A', 12, 103),
        ('1luc', 'A', 12, 158),
        #('1l3s', 'A', 12, 322),
        ('1usl', 'D', 12, 35),
        ('2fnu', 'B', 12, 179),
        ('1eu1', 'A', 12, 159),
        ('1ds1', 'A', 12, 28),
        ]

structures_16res = [
        ('1vzy', 'A', 16, 90), 
        ('1wdp', 'A', 16, 399),
        ('2z30', 'A', 16, 346),
        ('1izc', 'A', 16, 12), 
        ('1jq5', 'A', 16, 306),
        ('2iq7', 'D', 16, 110),
        ('1y0e', 'A', 16, 144),
        ('2rbd', 'B', 16, 72), 
        #('1kcz', 'B', 16, 70),
        ('2vjq', 'B', 16, 153),
        ('2v27', 'A', 16, 104),
        ('1or0', 'D', 16, 242),
        ('1vps', 'D', 16, 177),
        ('1h1n', 'A', 16, 276),
        ('1lvw', 'B', 16, 7), 
        ('2g8j', 'A', 16, 296),
        ('1h1n', 'B', 16, 201),
        ('1wdv', 'B', 16, 92),
        ('1vem', 'A', 16, 289),
        ('2vjq', 'A', 16, 389),
        ('2isa', 'C', 16, 94), 
        ('2aka', 'B', 16, 53),
        ('1ofw', 'B', 16, 17), 
        ('2cho', 'B', 16, 545),
        ('1t3t', 'A', 16, 617),
        ('1itw', 'A', 16, 704),
        ('1gte', 'D', 16, 318),
        ('1im5', 'A', 16, 95),
        ('1ddj', 'C', 16, 7),
        ]

structures_20res = [
        ('1w8o', 'A', 20, 151),
        ('2d81', 'A', 20, 301),
        ('1h4r', 'A', 20, 151), 
        ('1mj5', 'A', 20, 62),
        ('2q0z', 'X', 20, 237), 
        ('2afw', 'A', 20, 321), 
        ('1y7w', 'A', 20, 252), 
        ('2qp2', 'A', 20, 270),
        ('1m0z', 'A', 20, 179), 
        ('1vr5', 'B', 20, 359), 
        ('1oyc', 'A', 20, 289), 
        ('2isa', 'B', 20, 270),
        ('3bix', 'A', 20, 509), 
        ('1f2t', 'A', 20, 49), 
        ('1g2o', 'C', 20, 142), 
        ('1q6z', 'A', 20, 274),
        ('1ux6', 'A', 20, 890), 
        ('2cn3', 'A', 20, 395), 
        ('3bdi', 'A', 20, 62), 
        ('1w2l', 'A', 20, 21),
        ('2o7i', 'A', 20, 13), 
        ('1y7b', 'B', 20, 268), 
        ('1dhk', 'A', 20, 295), 
        ('1qhf', 'A', 20, 115),
        ('2ez9', 'A', 20, 544), 
        ('1elv', 'A', 20, 496), 
        ('1p4k', 'C', 20, 345), 
        ('1ne7', 'D', 20, 138),
        ('1y0p', 'A', 20, 48),
        ('3bi1', 'A', 20, 539),
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

def prepare_inputs(pdb_id, chain_to_keep, loop_size, loop_start):
    '''Prepare for the input files.'''
    download_pdb(pdb_id)

    # Extract the chains to keep

    raw_pose = rosetta.core.pose.Pose()
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(raw_pose, pdb_id + '.pdb')
   
    seqposes_to_keep = rosetta.utility.vector1_unsigned_long()
    for seqpos in range(1, raw_pose.size()):
        if raw_pose.pdb_info().chain(seqpos) == chain_to_keep:
            seqposes_to_keep.append(seqpos)
    
    rosetta.core.pose.pdbslice(pose, raw_pose, seqposes_to_keep)

    # Get the loop

    print pdb_id
    loop_seqpos_start = pose.pdb_info().pdb2pose(chain_to_keep, loop_start)
    loop = (loop_seqpos_start, loop_seqpos_start + loop_size - 1)

    with open(os.path.join('inputs', pdb_id + '.loop'), 'w') as f:
        f.write('LOOP {0} {1} {1} 0 1'.format(loop[0], loop[1]))

    # Generate the fasta file, this is important because the chains to be modeled are not labeled as chain A

    with open(os.path.join('inputs', pdb_id + '.fasta'), 'w') as f:
        f.write('>{0}|A\n'.format(pdb_id))
        f.write(pose.sequence())

    pose.dump_file(os.path.join('inputs', pdb_id + '.pdb'))
    pose.dump_file(os.path.join('reference', pdb_id + '.pdb'))

    os.remove(pdb_id + '.pdb')

def generate_dataset(set_name, structures):
    cwd = os.getcwd()
    
    if not os.path.exists(set_name):
        os.mkdir(set_name)
    os.chdir(set_name)

    if not os.path.exists('inputs'):
        os.mkdir('inputs')
    if not os.path.exists('reference'):
        os.mkdir('reference')
   
    for structure in structures:
        prepare_inputs(structure[0], structure[1], structure[2], structure[3])

    os.chdir(cwd)


if __name__ == '__main__':
    pyrosetta.init(options='-ignore_unrecognized_res')

    #generate_dataset('12res', structures_12res)
    generate_dataset('16res', structures_16res)
    generate_dataset('20res', structures_20res)
