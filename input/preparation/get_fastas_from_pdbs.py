#!/usr/bin/env python2

import os

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':

    pyrosetta.init()

    for f in os.listdir('.'):
        if f.endswith('.pdb'):
            pose = rosetta.core.pose.Pose()
            rosetta.core.import_pose.pose_from_file(pose, f)

            with open(f[:4] + '.fasta', 'w') as fout:
                fout.write('>' + f[:4] + '|A\n')
                sequence = pose.sequence()

                start = 0
                stop = start + 80

                while start < len(sequence):
                    fout.write(sequence[start:min(stop, len(sequence))] + '\n')
                    start = stop
                    stop = stop + 80

