#!/usr/bin/env python2

import gzip

from . import install
install.require_biopython()
import numpy as np
import Bio.PDB

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

def get_align_transformation_for_two_chains(chain1, chain2, residue_list1, residue_list2):
    '''Get the transformation that align chain1 to chain2 by 
    residues in two lists. Return the transformation matrix and vector
    '''
    assert(len(residue_list1) == len(residue_list2))
    
    def get_bb_coords_by_rosetta_ids(chain, residue_list):
        chain_list = [residue for residue in chain]
        valid_atoms = ['N', 'CA', 'C', 'O' ] #Only include the backbone atoms
        coords = []
        for residue in residue_list:
            for atom in chain_list[residue - 1]:
                if atom.get_name() in valid_atoms:
                    coords.append(atom.get_coord())

        return coords
    
    coords1 = get_bb_coords_by_rosetta_ids(chain1, residue_list1)
    coords2 = get_bb_coords_by_rosetta_ids(chain2, residue_list2)
    
    return get_superimpose_transformation(coords1, coords2)

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
            for atom in chain1_list[residue - 1]:
                if atom.get_name() in valid_atoms:
                    c1.append( atom.get_coord() )
            for atom in chain2_list[residue - 1]:
                if atom.get_name() in valid_atoms:
                    c2.append( atom.get_coord() )
        self.coord1 = np.array(c1)
        self.coord2 = np.array(c2)

    def rmsd(self, transformation=None):
        '''Calculate the RMSD between the coordinates. If a transformation
        defined as (M, t) is given, apply the transformation to coords1
        before calculating RMSD.
        '''
        if transformation is None:
            diff = self.coord1 - self.coord2
        else:
            M, t = transformation
            transformed_coord1 = np.array([np.dot(M, c) + t for c in self.coord1])
            diff = transformed_coord1 - self.coord2

        l = diff.shape[0]
        return np.sqrt( sum(sum( np.multiply(diff,diff) ))/l )


def calc_rmsd_from_file( file1, file2, residue_list, model1, model2, chain1_id=None, chain2_id=None, align_residues1=None, align_residues2=None):
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
    
    if align_residues1 is None:
        return rmsd_calculator.rmsd()
    else:
        M, t = get_align_transformation_for_two_chains(chain1, chain2, align_residues1, align_residues2)
        return rmsd_calculator.rmsd((M, t))
