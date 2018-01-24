#!/usr/bin/env python2

import gzip

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
        valid_atoms = ['N', 'CA', 'C', 'O' ] #Only include the backbone atoms
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


def calc_rmsd_from_file( file1, file2, residue_list, model1, model2, chain1_id=None, chain2_id=None):
    parser = Bio.PDB.PDBParser()
    
    def load_structure_file(sf):
        if sf.endswith('.pdb.gz'):
            with gzip.open(sf, 'r') as f:
                return parser.get_structure('', f)
        return parser.get_structure('', sf)

    # Load the structures
   
    structure1 = load_structure_file(file1)
    structure2 = load_structure_file(file2)

    # Get the chains

    if chain1_id is None:
        chain1 = [c for c in structure1[model1]][0]
    else:
        chain1 = structure1[model1][chain1_id]
    if chain2_id is None:
        chain2 = [c for c in structure2[model2]][0]
    else:
        chain2 = structure2[model2][chain2_id]

    rmsd_calculator = RMSDCalculator(chain1, chain2, residue_list)
    return rmsd_calculator.rmsd()
