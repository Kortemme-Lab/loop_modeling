#!/usr/bin/env python2

import numpy as np
import Bio.PDB

class RMSDCalculator():
    '''Calculate the RMSD of a given sequence between two structures''' 
    def __init__(self, chain1, chain2, residue_list):
        c1 = []
        c2 = []
        for residue in residue_list:
            for atom in chain1[residue]:
                if not atom.get_name().startswith('H'): #Only include heavy atoms
                    c1.append( atom.get_coord() )
            for atom in chain2[residue]:
                if not atom.get_name().startswith('H'): #Only include heavy atoms
                    c2.append( atom.get_coord() )
        self.coord1 = np.array(c1)
        self.coord2 = np.array(c2)
         
    def rmsd(self):
        diff = self.coord1 - self.coord2
        l = diff.shape[0]
        return np.sqrt( sum( sum(diff) )/l )

def calc_rmsd_from_file( file1, file2, model1, model2, chain1, chain2, residue_list ):
    parser = Bio.PDB.PDBParser()
    structure1 = parser.get_structure('', file1)
    structure2 = parser.get_structure('', file2)
    rmsd_calculator = RMSDCalculator(structure1[model1][chain1], structure2[model2][chain2], residue_list)
    return rmsd_calculator.rmsd()
