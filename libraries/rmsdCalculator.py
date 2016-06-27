#!/usr/bin/env python2

from . import install
install.require_biopython()
import numpy as np
import Bio.PDB

class RMSDCalculator():
    '''Calculate the RMSD of a given sequence between two structures''' 
    def __init__(self, chain1, chain2, residue_list):
        #Save residues of each chain into a list, such that residue indices become Rosetta numbers
        chain1_list = [residue for residue in chain1]
        chain2_list = [residue for residue in chain2]
        
        c1 = []
        c2 = []
        valid_atoms = ['N', 'CA', 'C' ] #Only include the backbone atoms
        for residue in residue_list:
            for atom in chain1_list[residue]:
                if atom.get_name() in valid_atoms:
                    c1.append( atom.get_coord() )
            for atom in chain2_list[residue]:
                if atom.get_name() in valid_atoms:
                    c2.append( atom.get_coord() )
        self.coord1 = np.array(c1)
        self.coord2 = np.array(c2)
         
    def rmsd(self):
        diff = self.coord1 - self.coord2
        l = diff.shape[0]
        return np.sqrt( sum(sum( np.multiply(diff,diff) ))/l )

def calc_rmsd_from_file( file1, file2, model1, model2, chain1, chain2, residue_list ):
    parser = Bio.PDB.PDBParser()
    structure1 = parser.get_structure('', file1)
    structure2 = parser.get_structure('', file2)
    rmsd_calculator = RMSDCalculator(structure1[model1][chain1], structure2[model2][chain2], residue_list)
    return rmsd_calculator.rmsd()
