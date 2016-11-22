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
from rmsdCalculator import calc_rmsd_from_file

class DataController:
    '''
    DataController controls the IO of data. 
    The lower level implementation could either use databases or dump data to files directly.
    '''
    def __init__(self, platform='disk'):
        self.platform = platform
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
    

    def calc_rmsd(self, loop_file, native_file, modeled_file):
        return self.implement.calc_rmsd(loop_file, native_file, modeled_file)
    

    def write_log(self, benchmark_id, protocol_id, stdout, stderr, job_id=0):
        return self.implement.write_log(benchmark_id, protocol_id, stdout, stderr, job_id) 
      

    def get_benchmark_list_by_name(self, database_name):
        return self.implement.get_benchmark_list_by_name(database_name)  

    def get_progress(self, database_name, benchmark_name):
        return self.implement.get_progress(database_name, benchmark_name)
    
    def get_unfinished_task_list(self, benchmark_id, num_all_tasks):
        return self.implement.get_unfinished_task_list(benchmark_id, num_all_tasks)

    def create_task_completion_list_file(self, benchmark_id, unfinished_task_list): 
        if self.platform == 'disk':
            return self.implement.create_task_completion_list_file(benchmark_id, unfinished_task_list)
    
    def read_task_completion_list(self, benchmark_id):
        if self.platform == 'disk':
            return self.implement.read_task_completion_list(benchmark_id)


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


    def get_benchmark_list_by_name(self, database_name):
        with database.connect(db_name = database_name) as session:
            return [r.name for r in session.execute('SELECT DISTINCT name from benchmarks ORDER BY benchmark_id DESC')]
       

    def get_progress(self, database_name, benchmark_name):
        # Create an entry in the benchmarks table.
        with database.connect(db_name = database_name) as session:
    
            messages = ['']
    
            # Use the latest benchmark's name if none was supplied
            if not benchmark_name:
                q = session.query(database.Benchmarks).order_by(database.Benchmarks.benchmark_id.desc())
                if q.count() == 0:
                    exit('There is no benchmark data in the database "{0}".'.format(database_name))
                benchmark_name = q.first().name
                messages.append('No benchmark was selected. Choosing the most recent benchmark: "{0}".\n'.format(benchmark_name))
    
            # Retrieve the set of benchmark runs associated the benchmark name
            q = session.query(database.Benchmarks).filter(database.Benchmarks.name == benchmark_name)
            if q.count() == 0:
                exit('There is no benchmark data in the database "{0}" for benchmark "{1}".'.format(database_name, benchmark_name))
            benchmark_runs = q.all()
    
            # Set nstruct to be the maximum value over the runs
            nstruct = max([b.nstruct for b in benchmark_runs])
    
            # Retrieve the set of PDB paths
            pdb_paths = set()
            for b in benchmark_runs:
                pdb_paths = pdb_paths.union(set([q.pdb_path for q in session.query(database.BenchmarkInputs).filter(database.BenchmarkInputs.benchmark_id == b.benchmark_id).all()]))
    
            # Retrieve the number of jobs with structures
            num_cases = len(pdb_paths)
            benchmark_size = nstruct * num_cases
            total_count = 0.0
            pdb_counts = dict.fromkeys(pdb_paths, 0)
            for r in session.execute('''
                    SELECT input_tag, COUNT(input_tag)
                    FROM structures
                    INNER JOIN batches ON structures.batch_id=batches.batch_id
                    INNER JOIN benchmark_protocols ON benchmark_protocols.protocol_id=batches.protocol_id
                    INNER JOIN benchmarks ON benchmarks.benchmark_id=benchmark_protocols.benchmark_id
                    WHERE benchmarks.name="{0}"
                    GROUP BY input_tag'''.format(benchmark_name)):
                total_count += min(r[1], nstruct) # do not count extra jobs
                assert(r[0] in pdb_paths)
                pdb_counts[r[0]] = r[1]
    
            num_failed = session.execute('''
                    SELECT COUNT(log_id) AS NumFailed
                    FROM tracer_logs
                    INNER JOIN benchmarks ON benchmarks.benchmark_id=tracer_logs.benchmark_id
                    WHERE stderr <> "" AND benchmarks.name="{0}"'''.format(benchmark_name))
            for r in num_failed: num_failed = r[0]; break
            progress = 100 * (total_count/benchmark_size)
    
            return dict(
                Messages = '\n'.join(messages),
                Progress = progress,
                StructureCount = num_cases,
                nstruct = nstruct,
                TotalCount = benchmark_size,
                CompletedCount = int(total_count),
                FailureCount = num_failed,
                CountPerStructure = pdb_counts,
            )
    
    
    def get_unfinished_task_list(self, benchmark_id, num_all_tasks):
        return None
        
    
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
        benchmark_define_dict['structures_path'] = os.path.join(self.data_path, id, 'structures')
        os.makedirs( os.path.join( self.data_path, id ) )
        os.makedirs( os.path.join( self.data_path, id, 'logs' ) )
        os.makedirs( os.path.join( self.data_path, id, 'structures' ) )
        
        # Dump the benchmark definition dictionary into a JSON file
        
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
    

    def calc_rmsd(self, loop_file, native_file, modeled_file):
        residue_list = []
        with open(loop_file, 'r') as f_loop:
            for line in f.readlines():
                s = line().split()
                for i in range(int(s[1]), int(s[2])+1):
                    residue_list.append(i) 
        self.rmsd = calc_rmsd_from_file(native_file, modeled_file, 0, 0, 'A', 'A', residue_list) 
    

    def write_log(self, benchmark_id, protocol_id, stdout, stderr, job_id):
        build_time_match = re.search(r"protocols.loop_modeling.LoopModeler: Build Time: ([\d]+) sec", stdout)
        build_time = int(build_time_match.groups()[0]) if build_time_match else 0
        
        centroid_time_match = re.search(r"protocols.loop_modeling.LoopModeler: Centroid Time: ([\d]+) sec", stdout)
        centroid_time = int(centroid_time_match.groups()[0]) if centroid_time_match else 0
        
        fullatom_time_match = re.search(r"protocols.loop_modeling.LoopModeler: Fullatom Time: ([\d]+) sec", stdout)
        fullatom_time = int(fullatom_time_match.groups()[0]) if fullatom_time_match else 0
        
        score_match = re.search(r"protocols.loop_modeling.LoopModeler: Total Score: ([-]?[\d\.]+)", stdout)
        if not score_match:
            score_match = re.search(r"protocols.loop_build.LoopBuildMover: total_energy ([-]?[\d\.]+)", stdout) #Legacy KIC

        score = float(score_match.groups()[0]) if score_match else None
        
        structure_match = re.search(r"-in:file:s ([\w\/\.]+)", stdout)
        structure = os.path.split(structure_match.groups()[0])[1].split('.')[0] if structure_match else None

        log_dict = { 'benchmark_id':benchmark_id,
                     'protocol_id':protocol_id,
                     'rmsd':self.rmsd,
                     'score':score,
                     'stderr':stderr,
                     'stdout':stdout,
                     'structure':structure,
                     'time_build':build_time,
                     'time_centroid':centroid_time,
                     'time_fullatom':fullatom_time}
        
        # Write the log dictionary into a JSON file and write the result atomically

        log_file = os.path.join( self.data_path, str(benchmark_id), 'logs', 'log'+str(job_id)+'.json') 
        with open( log_file, 'w' ) as f:
            f.write( json.dumps(log_dict, sort_keys=True, indent=4) )
        try:
            self.write_result(benchmark_id, job_id, structure, self.rmsd, score, build_time+centroid_time+fullatom_time)
        except:
            os.remove(log_file)
            raise

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
    

    def get_benchmark_list_by_name(self, database_name):
        benchmark_list = []
        for benchmark_id in os.listdir( self.data_path ):
            benchmark_define_dict = self.get_benchmark_define_dict( int(benchmark_id) )
            benchmark_list.append( (benchmark_id, benchmark_define_dict['name']) )
        return benchmark_list
       

    def get_progress(self, database_name, benchmark_name):
        progress_dict = {}
        if not benchmark_name:
            exit( "No benchmark was selected!" )

        #Get the benchmark ids coresponding to the benchmark name

        all_benchmark_list = self.get_benchmark_list_by_name(database_name)
        benchmark_list = [ int(entry[0]) for entry in all_benchmark_list if entry[1] == benchmark_name ]
        most_recent_id = max( benchmark_list )
        progress_dict['MostRecentID'] = most_recent_id
        benchmark_define_dict = self.get_benchmark_define_dict(most_recent_id)

        #Get the Messages

        progress_dict['Messages'] = 'Reporting the most recent benchmark with name {0}, id={1}'.format(benchmark_name, most_recent_id)

        #Get nstruct

        progress_dict['nstruct'] = benchmark_define_dict['nstruct']
        pdb_pathes = set( [ i.pdb_path for i in benchmark_define_dict['input_pdbs'] ] )
        pdb_name_path_dict = {}
        for p in pdb_pathes:
            pdb_name_path_dict[ os.path.split(p)[1].split('.')[0] ] = p
        progress_dict['StructureCount'] = len(pdb_pathes)
        progress_dict['TotalCount'] = progress_dict['nstruct'] * progress_dict['StructureCount']

        #Count completed jobs

        progress_dict['CountPerStructure'] = dict.fromkeys(pdb_pathes, 0)
        progress_dict['CompletedCount'] = 0
        with open( os.path.join(self.data_path, str(most_recent_id), benchmark_name+'.results'), 'r' ) as f_result:
            for line in f_result.readlines():
                if line.startswith('#'): continue
                s = line.split()
                progress_dict['CountPerStructure'][ pdb_name_path_dict[s[0]] ] += 1
                progress_dict['CompletedCount'] += 1 
        progress_dict['Progress'] = 100 * progress_dict['CompletedCount']/progress_dict['TotalCount']

        #Count failed jobs

        progress_dict['FailureCount'] = 0
        for log in os.listdir( os.path.join(self.data_path, str(most_recent_id), 'logs') ):
            with open( os.path.join(self.data_path, str(most_recent_id), 'logs', log) ) as f:
                log_dict = json.loads( f.read() )  
                if log_dict['stderr'] != '':
                    progress_dict['FailureCount'] += 1
        return progress_dict


    def get_unfinished_task_list(self, benchmark_id, num_all_tasks):
        log_list = os.listdir( os.path.join(self.data_path, str(benchmark_id), 'logs') )
        finished_task_list = [ int(l.split('.')[0][3:]) for l in log_list ] 
        unfinished_task_list = list( set(range(0,num_all_tasks)).difference( set(finished_task_list) ) )
        return unfinished_task_list

    def create_task_completion_list_file(self, benchmark_id, unfinished_task_list): 
        '''Create a file that contains the list of unfinished tasks'''
        with open(os.path.join(self.data_path, str(benchmark_id), 'task_completion_list'), 'w') as f:
            for task_num in unfinished_task_list:
                f.write(str(task_num)+'\n')

    def read_task_completion_list(self, benchmark_id):
        '''Read the task completion list file into a list'''
        with open(os.path.join(self.data_path, str(benchmark_id), 'task_completion_list'), 'r') as f:
            return [ int(task_num) for task_num in f.readlines() ]


class DiskBenchmarkInput:
    ''' For the backward compatibility to the database.BenchmarkInput '''
    def __init__ (self, pdb_path):
        self.pdb_path = pdb_path
