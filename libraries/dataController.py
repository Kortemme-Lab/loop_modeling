#!/usr/bin/env python2

import getpass
import json
import database

class DataController:
    '''
    DataController controls the IO of data. 
    The lower level implementation could either use databases or dump data to files directly.
    '''
    def __init__(self, platform='disk'):
        if 'disk' == platform:
            self.implement = DiskDataController()
        elif 'database' == platform:
            self.implement = DatabaseDataController()
        else:
            raise ValueError('Invalid dataController platform %s' % platform )
    
    def create_benchmark(self, benchmark_define_dict, pdbs):
        benchmark_id = self.implement.create_benchmark(benchmark_define_dict, pdbs)
        print( "Your benchmark \"{0}\" (id={1}) has been created".format(
               benchmark_define_dict['name'], benchmark_id) )
        return benchmark_id
        
    def get_benchmark_define_dict(self, benchmark_id):
        return self.implement.get_benchmark_define_dict(benchmark_id)
        
    def get_benchmark_variables(self, benchmark_id):
        return self.implement.get_benchmark_variables(benchmark_id)
    
    def write_log(self, benchmark_id, protocol_id, stdout, stderr):
        return self.implement.write_log(benchmark_id, protocol_id) 
      
class DatabaseDataController:
    '''
    Implementation of DataController on databases.
    '''
    def __init__(self):
        # Test the database connection
        try: database.test_connect()
        except RuntimeError, error: 
            print error
            sys.exit(1)

    def create_benchmark(self, benchmark_define_dict, pdbs):
        d = benchmark_define_dict
        # Create an entry in the benchmarks table.
        with database.connect() as session:
            benchmark = database.Benchmarks(
                    d['name'], d['script'], d['nstruct'],
                    user=d['user'], desc=d['desc'],
                    vars=d['vars'], flags=d['flags'], fragments=d['fragments'],
                    git_commit=d['git_commit'], git_diff=d['git_diff'],
                    fast=d['fast'], non_random=d['non_random'],
            )

            for pdb in pdbs:
                benchmark_input = database.BenchmarkInputs(pdb)
                benchmark.input_pdbs.append(benchmark_input)

            session.add(benchmark); session.flush()
            return str(benchmark.id)

    def get_benchmark_define_dict(self, benchmark_id):
        benchmark_define_dict = {}
        with database.connect() as session:
            benchmark = session.query(database.Benchmarks).get(benchmark_id)
            benchmark_define_dict['id'] = benchmark.id
            benchmark_define_dict['name'] = benchmark.name 
            benchmark_define_dict['script'] = benchmark.rosetta_script
            benchmark_define_dict['nstruct'] = benchmark.nstruct
            benchmark_define_dict['user'] = benchmark.user
            benchmark_define_dict['desc'] = benchmark.description
            benchmark_define_dict['vars'] = benchmark.rosetta_script_vars
            benchmark_define_dict['flags'] = benchmark.rosetta_flags 
            benchmark_define_dict['fragments'] = benchmark.rosetta_fragments 
            benchmark_define_dict['git_commit'] = benchmark.git_commit
            benchmark_define_dict['git_diff'] = benchmark.git_diff
            benchmark_define_dict['fast'] = benchmark.fast
            benchmark_define_dict['non_random'] = benchmark.non_random
            benchmark_define_dict['input_pdbs'] = benchmark.input_pdbs
        return benchmark_define_dict

    def get_benchmark_variables(self, benchmark_id):
        benchmark_records = [r for r in session.query(database.Benchmarks).filter(database.Benchmarks.name == benchmark_id)]
        benchmark_variables = dict(
            rosetta_script = set([r.rosetta_script for r in benchmark_records]),
            rosetta_script_vars = [json.loads(r.rosetta_script_vars) for r in benchmark_records],
            rosetta_flags = set([r.rosetta_flags for r in benchmark_records]),
            rosetta_fragments = set([r.rosetta_fragments for r in benchmark_records]),
            fast = set([r.fast for r in benchmark_records]),
            non_random = set([r.non_random for r in benchmark_records]),
        )
        return benchmark_variables

    def write_log(self, benchmark_id, protocol_id, stdout, stderr):
        with database.connect() as session:
            if protocol_id is not None:
                benchmark_map = database.BenchmarkProtocols(benchmark_id, protocol_id)
                session.add(benchmark_map)
                session.commit()  # Make sure the protocol mapping is saved even if 
                                  # something else messes up this transaction later on.
            log_row = database.TracerLogs(benchmark_id, protocol_id, stdout, stderr)
            session.add(log_row)

    
class DiskDataController:
    '''
    Implementation of dataController on disk.
    '''
    pass
