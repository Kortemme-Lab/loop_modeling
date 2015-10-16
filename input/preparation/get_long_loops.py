import os

from tools.bio.rcsb import retrieve_pdb
from tools.fs.fsio import write_file

tbl = '''
1o97	1.6A	156-169
1ock	1.8A	209-222
1r6x	1.40A	72-85
1ra0	1.12A	361-375
1rdq	1.26A	272-285
1rv9	1.53A	225-238
1s95	1.60A	477-491
1xu1	1.90A	221-234
1zhx	1.50A	392-406
2aeb	1.29A	156-170
2b0t	1.75A	701-715
2bwr	1.5A	269-282
2c0h	1.6A	40-53
2cjp	1.95A	58-72
2hdw	2.00A	131-147
2o2k	1.60A	1220-1234
2o2k	1.60A	1221-1234
2puh	1.82A	70-85
2v3v	1.99A	382-396
3a3p	1.90A	286-300
3a64	1.60A	350-364
3b40	2.00A	389-402
3bb7	1.50A	231-245
3cnq	1.71A	50-63
3css	1.70A	95-109
3ea1	1.75A	136-150
3f1l	0.95A	98-112
3h2g	1.86A	124-140
'''

cases = [l.strip().split() for l in tbl.split('\n') if l.strip()]
for c in cases:
    pdb_store = '../structures/14_17_res/rcsb/reference'
    pdb_id = c[0]
    pdb_path = os.path.join(pdb_store, '{0}.pdb'.format(pdb_id))
    write_file(pdb_path, retrieve_pdb(pdb_id))

