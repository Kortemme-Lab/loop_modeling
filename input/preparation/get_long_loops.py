import os
import pprint
import json

from tools import colortext
from tools.bio.rcsb import retrieve_pdb
from tools.bio.pdb import PDB
from tools.fs.fsio import read_file, write_file

# This list was compiled from the supplementals in Stein & Kortemme (doi:10.1371/journal.pone.0063090) and referenced with
# Zhao et al. (doi:10.1002/prot.23129) to get the chain information
loop_definitions = '''
1o97	1.6A	D:156-D:169
1ock	1.8A	A:209-A:222
1r6x	1.40A	A:72-A:85
1ra0	1.12A	A:361-A:375
1rdq	1.26A	E:272-E:285
1rv9	1.53A	A:225-A:238
1s95	1.60A	A:477-A:491
1xu1	1.90A	A:221-A:234
1zhx	1.50A	A:392-A:406
2aeb	1.29A	B:156-B:170
2b0t	1.75A	A:701-A:715
2bwr	1.5A	A:269-A:282
2c0h	1.6A	A:40-A:53
2cjp	1.95A	A:58-A:72
2hdw	2.00A	A:131-A:147
#2o2k	1.60A	A:1220-A:1234 # Note: this is too similar to the next entry so we discard it
2o2k	1.60A	A:1221-A:1234
2puh	1.82A	A:70-A:85
2v3v	1.99A	A:382-A:396
3a3p	1.90A	A:286-A:300
3a64	1.60A	A:350-A:364
3b40	2.00A	A:389-A:402
3bb7	1.50A	A:231-A:245
3cnq	1.71A	S:50-S:63
3css	1.70A	A:95-A:109
3ea1	1.75A	A:136-A:150
3f1l	0.95A	A:98-A:112  # Note: this was A:99-A:113 in Zhao et al.
3h2g	1.86A	A:124-A:140
'''

def setup_pdb_files_and_loop_definitions():
    print('Downloading files')
    bcases = [l.strip().split() for l in loop_definitions.split('\n') if l.strip() and l[0] != '#']
    n = 0
    loop_definitions = {}
    for bcase in bcases:
        n += 1
        pdb_id = bcase[0]
        colortext.message('\nSetting up files for case #{0}/{1}: {2}'.format(n, len(bcases), pdb_id))
        pdb_store = os.path.join('..', 'structures', '14_17_res', 'rcsb', 'reference')
        pdb_path = os.path.join(pdb_store, '{0}.pdb'.format(pdb_id))
        if not os.path.exists(pdb_path):
            write_file(pdb_path, retrieve_pdb(pdb_id))

        p = PDB.from_filepath(pdb_path)
        assert(p.seqres_sequences.keys() == p.atom_sequences.keys())
        if len(p.seqres_sequences) > 1:
            num_unique_sequences = len(set(map(str, p.seqres_sequences.values())))
            if num_unique_sequences == 1:
                colortext.pcyan('homo {0}-mer'.format(len(p.seqres_sequences)))
            elif num_unique_sequences < len(p.seqres_sequences):
                colortext.plightpurple('some duplicated sequences')
            else:
                colortext.warning('possible heteromer')
        else:
            assert(len(p.seqres_sequences) == 1)
            colortext.pcyan('monomer')

        for c, seq in sorted(p.seqres_sequences.iteritems()):
            print('{0}: {1}'.format(c, seq))

        assert(bcase[1].endswith('A'))
        assert(float(bcase[1][:-1]) == p.get_resolution())

        assert(len(bcase[2].split('-')) == 2)
        loop_start = bcase[2].split('-')[0]
        loop_end = bcase[2].split('-')[1]
        assert(loop_start[0] == loop_end[0] and loop_start[1] == ':' and loop_end[1] == ':')
        chain_id = loop_start[0]
        start_residue = str(int(loop_start[2:]))
        end_residue = str(int(loop_end[2:]))

        loop_sequence = ''
        in_loop = False
        for respair in p.atom_sequences[chain_id]:
            if str(respair[1]) == PDB.ChainResidueID2String(chain_id, start_residue):
                in_loop = True
            if in_loop:
                loop_sequence += respair[1].ResidueAA
            if str(respair[1]) == PDB.ChainResidueID2String(chain_id, end_residue):
                break
        assert(14 <= len(loop_sequence) <= 17)

        #import sys
        #sys.exit(0)
        assert(pdb_id not in loop_definitions)
        loop_definitions[pdb_id] = dict(
            chainID = chain_id,
            DOI = '10.1002/prot.23129',
            EndResidueID = PDB.ResidueID2String(end_residue),
            PassedFilter = True,
            Sequence = loop_sequence,
            StartResidueID = PDB.ResidueID2String(start_residue),
        )

    loop_definitions_path = os.path.join('..', 'structures', '14_17_res', 'loop_definitions.json')
    if not os.path.exists(loop_definitions_path):
        write_file(loop_definitions_path, json.dumps(loop_definitions, indent = 4, sort_keys = True))


def create_loop_json_files():
    loop_definitions = json.loads(read_file(os.path.join('..', 'structures', '14_17_res', 'loop_definitions.json')))

    for pdb_id, details in sorted(loop_definitions.iteritems()):

        # The residues happen to not have insertion codes
        EndResidueID = int(details['EndResidueID'])
        StartResidueID = int(details['EndResidueID'])

        loop_json = {
            "LoopSet": [
                {
                    "cut": None,
                    "extras": {
                        "extend": None,
                        "skip_rate": None
                    },
                    "start": {
                        "chainID": details['chainID'],
                        "iCode": " ",
                        "resSeq": StartResidueID
                    },
                    "stop": {
                        "chainID": details['chainID'],
                        "iCode": " ",
                        "resSeq": EndResidueID
                    }
                }
            ]
        }
        loop_json_filepath = os.path.join('..', 'structures', '14_17_res', 'rosetta', 'pruned', '{0}.loop.json'.format(pdb_id))
        if not os.path.exists(loop_json_filepath):
            write_file(loop_json_filepath, json.dumps(loop_json, sort_keys = True, indent = 4))

