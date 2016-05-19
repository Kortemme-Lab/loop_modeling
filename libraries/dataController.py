#!/usr/bin/env python2

import getpass
import json
import database
import os
import sys
import re
import fcntl
import time
import errno

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
    
    def write_log(self, benchmark_id, protocol_id, stdout, stderr, job_id=0):
        return self.implement.write_log(benchmark_id, protocol_id, stdout, stderr, job_id) 
      
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
                    vars=json.dumps(d['vars']), flags=d['flags'], fragments=d['fragments'],
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

    def write_log(self, benchmark_id, protocol_id, stdout, stderr, job_id):
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
    def __init__(self):
        self.data_path = './data' # Save the data into the './data' directory. Maybe make it user defined latter.
        if not os.path.exists(self.data_path):
            try: os.makedirs(self.data_path)
            except RuntimeError, error: 
                print error
                sys.exit(1)

    def create_benchmark(self, benchmark_define_dict, pdbs):
        id = self.get_new_benchmark_id()
        benchmark_define_dict['input_pdbs'] = pdbs
        benchmark_define_dict['id'] = id
        os.makedirs( os.path.join( self.data_path, id ) )
        os.makedirs( os.path.join( self.data_path, id, 'logs' ) )
        with open( os.path.join( self.data_path, id, 'benchmark_define.json' ), 'w' ) as f:
            f.write(json.dumps(benchmark_define_dict, sort_keys=True, indent=4) )
        with open( os.path.join( self.data_path, id, benchmark_define_dict['name']+'.results' ), 'w' ) as f:
            f.write('#PDB    Model   Loop_rmsd   Total_energy    Runtime\n')
        return id
    
    def get_new_benchmark_id(self):
        current_ids = os.listdir(self.data_path)
        return '0' if 0==len(current_ids) else str(1 + max( [ int(x) for x in current_ids] )) 

    def get_benchmark_define_dict(self, benchmark_id):
        benchmark_define_dict = None
        with open( os.path.join( self.data_path, str(benchmark_id), 'benchmark_define.json' ), 'r' ) as f:
            benchmark_define_dict = json.loads( f.read() )  
        benchmark_define_dict['input_pdbs'] = [DiskBenchmarkInput(pdb_path) for pdb_path in benchmark_define_dict['input_pdbs']]
        return benchmark_define_dict

    def get_benchmark_variables(self, benchmark_id, job_id):
        # TODO: implement
        return None
    
    def write_log(self, benchmark_id, protocol_id, stdout, stderr, job_id):
        build_time_match = re.search(r"protocols.loop_modeling.LoopModeler: Build Time: ([\d]+) sec", stdout)
        build_time = int(build_time_match.groups()[0]) if build_time_match else None
        centroid_time_match = re.search(r"protocols.loop_modeling.LoopModeler: Centroid Time: ([\d]+) sec", stdout)
        centroid_time = int(centroid_time_match.groups()[0]) if centroid_time_match else None
        fullatom_time_match = re.search(r"protocols.loop_modeling.LoopModeler: Fullatom Time: ([\d]+) sec", stdout)
        fullatom_time = int(fullatom_time_match.groups()[0]) if fullatom_time_match else None
        rmsd_match = re.search(r"protocols.loop_modeling.LoopModeler: Loop Backbone RMSD: ([\d\.]+)", stdout)
        rmsd = float(rmsd_match.groups()[0]) if rmsd_match else None
        score_match = re.search(r"protocols.loop_modeling.LoopModeler: Total Score: ([-]?[\d\.]+)", stdout)
        score = float(score_match.groups()[0]) if score_match else None
        structure_match = re.search(r"-in:file:s ([\w\/\.]+)", stdout)
        structure = os.path.split(structure_match.groups()[0])[1].split('.')[0] if structure_match else None

        log_dict = { 'benchmark_id':benchmark_id,
                     'protocol_id':protocol_id,
                     'rmsd':rmsd,
                     'score':score,
                     'stderr':stderr,
                     'stdout':stdout,
                     'structure':structure,
                     'time_build':build_time,
                     'time_centroid':centroid_time,
                     'time_fullatom':fullatom_time}
        with open( os.path.join( self.data_path, str(benchmark_id), 'logs', 'log'+str(job_id)+'.json'), 'w' ) as f:
            f.write( json.dumps(log_dict, sort_keys=True, indent=4) )
        self.write_result(benchmark_id, job_id, structure, rmsd, score, build_time+centroid_time+fullatom_time)

    def write_result(self, benchmark_id, job_id, structure, rmsd, score, runtime):
        benchmark_define_dict = self.get_benchmark_define_dict( benchmark_id )
        model_id = job_id // len(benchmark_define_dict['input_pdbs'])
        #Accquire the lock that protects the result file
        lock_fd = os.open( os.path.join( self.data_path, str(benchmark_id), 'lock' ), os.O_CREAT)
        while True:
            try:
                fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except IOError as e:
                if e.errno != errno.EAGAIN:
                    raise
                else:
                    time.sleep(0.1)
        try: 
            with open( os.path.join( self.data_path, str(benchmark_id), benchmark_define_dict['name']+'.results' ), 'a' ) as f:
                f.write('%s\t%d\t%f\t%f\t%d\n' % (structure, model_id, rmsd, score, runtime) )
        except:
            fcntl.flock(lock_fd, fcntl.LOCK_UN)
            raise
        #Release the lock
        fcntl.flock(lock_fd, fcntl.LOCK_UN)
        os.close(lock_fd) 

class DiskBenchmarkInput:
    ''' For the backward compatibility to the database.BenchmarkInput '''
    def __init__ (self, pdb_path):
        self.pdb_path = pdb_path
