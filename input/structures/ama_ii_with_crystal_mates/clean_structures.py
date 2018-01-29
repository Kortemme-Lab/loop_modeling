#!/usr/bin/env python2

import os

import pyrosetta
from pyrosetta import rosetta

if __name__ == '__main__':
    pyrosetta.init(options='-ignore_unrecognized_res')

    for f in os.listdir('.'):
        if f.endswith('.pdb'):

            pose = rosetta.core.pose.Pose()
            rosetta.core.import_pose.pose_from_file(pose, f)
            
            pose.dump_file(f)
