
def get_total_score_from_pdb(pdb_file):
    '''Get the rosetta total score from a pdb file'''
    with open(pdb_file, 'r') as f:
        for line in f.readlines():
            if line.startswith('total_score'):
                return float(line.split()[-1])
    return None


if __name__ == '__main__':
    import sys
    
    print(get_total_score_from_pdb(sys.argv[1]))
