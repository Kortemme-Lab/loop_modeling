#!/usr/bin/env python2
# encoding: utf-8

# The MIT License (MIT)
#
# Copyright (c) 2015 Shane O'Connor
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""\
remove_native_information.py
Functions to help remove particular side-chains or atoms from a PDB structure. This file does not typically need to be run;
it is included to document how the input structures were generated. There are two functions in this script:
 - create_pruned_structures removes the loop residues and surrounding sidechains (10A) for each dataset case from both the
   RCSB structures and the structures minimized using the Rosetta force field. The outputs are a PDB structure and a FASTA
   file of the removed loop sequence. We do this to remove native bias from the dataset;
 - add_missing_residues adds the residues in the FASTA file back into the Rosetta-minimized structures in a non-native
   conformation (backbone atoms are placed almost equidistantly in a straight line between the buttressing residues) and
   creates a Rosetta .loop file in Rosetta numbering (residues are numbered sequentially using integers starting from 1).
   This processing allows Rosetta to work with the structures but without adding back any native bias.
This script depends on the Kortemme Lab tools repository. This will be made available at https://github.com/Kortemme-Lab/tools
in the future.

Setup:
The tools repository will be added as a Git submodule but can be set up as below. All commands should be run in the
directory containing this script:

i) git clone https://github.com/Kortemme-Lab/tools.git
or
ii) (For Linux, Mac OS X, Cygwin): Clone the repository elsewhere and:
    ln -s /path/to/tools tools

Usage:
    remove_native_information <output_directory>

Arguments:
    <output_directory>
        The name of the directory where the output structures will be created.

Created by Shane O'Connor 2015.
"""

import re
import sys
import os
import types
import string
import types
import pprint
import math
import numpy
import glob
import json
import traceback

if __name__ == '__main__':
    sys.path.insert(0, os.path.join('..', '..'))

from tools import colortext
import tools.bio.rcsb
from tools.bio.basics import Residue, PDBResidue, Sequence, SequenceMap, residue_type_3to1_map, protonated_residue_type_3to1_map, non_canonical_amino_acids, protonated_residues_types_3, residue_types_3, Mutation, ChainMutation, SimpleMutation
from tools.bio.basics import dna_nucleotides, rna_nucleotides, dna_nucleotides_3to1_map, dna_nucleotides_2to1_map, non_canonical_dna, non_canonical_rna, all_recognized_dna, all_recognized_rna
from tools.bio.bonsai import Bonsai, PDBSection
from tools.bio.fasta import FASTA
from tools.bio.pdb import PDB
from tools.bio.spackle import Spackler
from tools.rosetta.map_pdb_residues import get_pdb_contents_to_pose_residue_map
from tools.fs.fsio import read_file, write_file
from tools.pymath.stats import get_mean_and_standard_deviation
from tools.pymath.cartesian import spatialhash
from tools.rosetta.map_pdb_residues import get_pdb_contents_to_pose_residue_map
from tools.general.strutil import remove_trailing_line_whitespace as normalize_pdb_file
from tools.rosetta.input_files import LoopsFile


def prepare_structures(file_filter, output_directory, loop_definitions, require_filter = True, create_partial_structures = False):
    search_radius = 10.0

    if not(os.path.exists(output_directory)):
        os.mkdir(output_directory)

    # Iterate through the dataset cases
    for pdb_file in sorted(glob.glob(file_filter)):
        pdb_prefix = os.path.splitext(os.path.split(pdb_file)[1])[0].lower()

        # Read the benchmark loop definition
        if not loop_definitions.get(pdb_prefix):
            raise Exception('The loop definition for {0} is missing.'.format(pdb_prefix))
        loop_definition = loop_definitions[pdb_prefix]
        loops = [PDBSection(loop_definition['chainID'], loop_definition['StartResidueID'], loop_definition['EndResidueID'], Sequence = loop_definition['Sequence'])]

        # Only process files that passed the benchmark criteria
        if require_filter and not loop_definition['PassedFilter']:
            continue

        # Remove the loops and surrounding sidechain atoms from the structure
        b = Bonsai(read_file(pdb_file))
        bonsai, cutting, PSE_file, PSE_script, FASTA_file = b.prune_loop_for_kic(loops, search_radius, expected_loop_length = 12, generate_pymol_session = True)

        # Create a PyMOL session file for visual inspection
        write_file(os.path.join(output_directory, '{0}.pse'.format(pdb_prefix)), PSE_file)

        # Create the new PDB file with the loop and surrounding sidechains removed
        write_file(os.path.join(output_directory, '{0}.pdb'.format(pdb_prefix)), bonsai)
        if create_partial_structures:
            write_file(os.path.join(output_directory, '{0}_missing_loop_and_surrounding_sidechains.pdb'.format(pdb_prefix)), bonsai)
            write_file(os.path.join(output_directory, '{0}_loop_and_surrounding_sidechains.pdb'.format(pdb_prefix)), cutting)

        # Create the FASTA file containing the loop sequence. This will be used along with the loop_definitions.json file
        # to add the residues back into the Rosetta structure
        write_file(os.path.join(output_directory, '{0}.fasta'.format(pdb_prefix)), FASTA_file)

        sys.stdout.write('.')
        sys.stdout.flush()

    print('')


def create_pruned_structures(output_directory):
    loop_definitions = json.loads(read_file('../structures/loop_definitions.json'))
    prepare_structures('../structures/12_res/rcsb/original/*.pdb', os.path.join(output_directory, '12_res_rcsb'), loop_definitions)
    prepare_structures('../structures/12_res/rosetta/preminimized/*.pdb', os.path.join(output_directory, '12_res_rosetta'), loop_definitions)


def add_missing_residues(output_directory):
    '''Add the removed loop residue backbone atoms back into the Rosetta preminimized structures but in a clearly non-native manner.'''

    output_directory = os.path.join(output_directory, 'rosetta')
    if not(os.path.exists(output_directory)):
        os.mkdir(output_directory)

    # Determine path to RosettaScripts
    rosetta_scripts_binary = None
    release_binaries = []
    other_binaries = []
    from libraries import settings
    try:
        settings.load(arguments, allow_failure = True)
    except Exception, e:
        raise colortext.Exception('An exception occurred reading the settings.conf file: {0}.'.format(str(e)))
    rosetta_binary_path = os.path.join(settings.rosetta, 'source', 'bin')
    for rosetta_scripts_binary in glob.glob(os.path.join(rosetta_binary_path, 'rosetta_scripts*')):
        if 'release' in rosetta_scripts_binary:
            release_binaries.append((len(rosetta_scripts_binary), rosetta_scripts_binary)) # shorter filenames usually indicate less specialized Rosetta binaries
        else:
            other_binaries.append((len(rosetta_scripts_binary), rosetta_scripts_binary))
    if release_binaries:
        rosetta_scripts_binary = sorted(release_binaries)[0]
    elif other_binaries:
        rosetta_scripts_binary = sorted(other_binaries)[0]
    if not rosetta_scripts_binary:
        raise colortext.Exception('No RosettaScripts binary could be located in {0}.'.format(rosetta_binary_path))
    rosetta_scripts_binary = rosetta_scripts_binary[1]

    # Iterate through the dataset cases
    file_filter = '../structures/rosetta/pruned/*.pdb'
    for pdb_file in sorted(glob.glob(file_filter)):
        pdb_prefix = os.path.splitext(os.path.split(pdb_file)[1])[0].lower()
        file_prefix = os.path.splitext(pdb_file)[0]
        fasta_file = file_prefix + '.fasta'
        loop_file = file_prefix + '.loop.json'
        assert(os.path.exists(fasta_file))
        assert(os.path.exists(loop_file))

        # Convert the FASTA headers back into PDB residue IDs
        fasta_contents = read_file(fasta_file)
        headers = [l for l in fasta_contents.split('\n') if l.startswith('>')]
        assert(len(headers) == 1)
        header = headers[0]
        pdb_residue_ids = [PDB.ChainResidueID2String(l[0], l[1:]) for l in header[header.find('Residues ') + 9:].split(';')]

        # Add the missing atoms atoms back into the PDB file
        spackler = Spackler.from_filepath(pdb_file)
        new_pdb_content = spackler.add_backbone_atoms_linearly_from_loop_filepaths(loop_file, fasta_file, pdb_residue_ids)
        write_file(os.path.join(output_directory, '{0}.pdb'.format(pdb_prefix)), new_pdb_content)

        # Create a Rosetta .loop file
        loop_set = json.loads(read_file(loop_file)).get('LoopSet')
        assert(len(loop_set) == 1)
        start_res = '{chainID}{resSeq:>4d}{iCode}'.format(**loop_set[0]['start'])
        end_res = '{chainID}{resSeq:>4d}{iCode}'.format(**loop_set[0]['stop'])

        success, result = get_pdb_contents_to_pose_residue_map(new_pdb_content, rosetta_scripts_binary, None, pdb_id = None, extra_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res')

        if not success:
            colortext.error('Failed on {0}.'.format(pdb_prefix))
            raise colortext.Exception('\n'.join(result))
        else:
            if not start_res in result:
                raise colortext.Exception('Could not find the starting residue in the PDB -> Rosetta residue mapping.')
            elif not end_res in result:
                raise colortext.Exception('Could not find the starting residue in the PDB -> Rosetta residue mapping.')
            start_rosetta_res = result[start_res]['pose_residue_id']
            end_rosetta_res = result[end_res]['pose_residue_id']
            if not end_rosetta_res > start_rosetta_res:
                raise colortext.Exception('The end residue have a higher index number than the starting residue.')
            loop_file_content = 'LOOP {0} {1}\n'.format(start_rosetta_res, end_rosetta_res)

            # Create the new PDB file with the loop residue backbone atoms added back in
            write_file(os.path.join(output_directory, '{0}.pdb'.format(pdb_prefix)), new_pdb_content)

            # Create a loop file in Rosetta numbering for protocols which require that format (the loopmodel code currently
            # requires this at the time of writing)
            write_file(os.path.join(output_directory, '{0}.loop'.format(pdb_prefix)), loop_file_content)

            # Remove this block after removing the old .loop files
            lfc = loop_file_content.strip().split()
            ofc = read_file(os.path.join('..', 'structures', 'rosetta', pdb_prefix + '.loop')).strip().split()
            assert(lfc[1] == ofc[0] and lfc[2] == ofc[1])

        sys.stdout.write('.')
        sys.stdout.flush()

    print('')


if __name__ == '__main__':
    from libraries import docopt
    arguments = docopt.docopt(__doc__)
    output_directory = arguments['<output_directory>']
    e, trc = '', ''

    if True:
        # Disable this code by default
        try:
            os.mkdir(output_directory)
        except Exception, e: trc = traceback.format_exc()
        if not os.path.exists(output_directory):
            colortext.error('Error: Could not create the output directory.')
            if e: colortext.error(str(e))
            colortext.warning(trc)

        #create_pruned_structures(output_directory)
        add_missing_residues(output_directory)


