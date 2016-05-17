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


class DiskDataController:
    '''
    Implementation of dataController on disk.
    '''
    pass
