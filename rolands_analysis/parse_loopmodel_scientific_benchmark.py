#!/usr/bin/python
#This file was developed and written by Roland A. Pache, Copyright (C) 2011, 2014.

#import libraries
import re
import os
import time
import sys
from subprocess import *
#from optparse import OptionParser
#import Parameters
import Statistics
import functions_lib
import datetime
#import Proteins


print '\nParses the loopmodel scientific benchmark.\n'


#constants
input_time_format="%a %b %d %H:%M:%S %Z %Y"
output_time_format="%s"

#patterns

#define classes
class Model:
    def __init__(self):
        self.id=None
        self.loop_rms=float('inf')
        self.total_energy=float('inf')
        self.runtime=None

#define functions
def getEnergyStats(models):
    energies=[]
    for model in models:
        if model.id!=None:
            energies.append(model.total_energy)
    #--
    return (min(energies),functions_lib.median(energies),max(energies))

def getRmsdStats(models):
    rmsds=[]
    for model in models:
        if model.id!=None:
            rmsds.append(model.loop_rms)
    #--
    return (min(rmsds),functions_lib.median(rmsds),max(rmsds))

            
start_time=time.time()

#parse input parameters
if len(sys.argv)!=2:
    print
    print 'Usage: ./parse_loopmodel_scientific_benchmark.py PARAMETER_FILE'
    print
    sys.exit()
#-
parameter_file=sys.argv[1]
parameters=functions_lib.parseParameterFile(parameter_file)
main_dir=parameters['global_main_dir']
num_models_offset=int(parameters['loopmodel_num_models_offset'])
num_models_per_PDB=int(parameters['loopmodel_num_models_per_PDB'])
run_identifier=parameters['loopmodel_run_identifier']
indir=main_dir+parameters['loopmodel_outdir']+run_identifier+'/'
outfile_name=main_dir+parameters['loopmodel_outdir']+run_identifier+'.results'

#prepare outfile
outfile=open(outfile_name,'w')
outfile.write('#PDB\tModel\tLoop_rmsd\tTotal_energy\tRuntime\n')

#parse benchmark results
print indir
sorted_indir_contents=sorted(os.listdir(indir))
best_models=[]
closest_models=[]
runtimes=[]
for pdb in sorted_indir_contents:
    print
    print pdb
    #parse cluster outfiles and collect model stats
    models=[]
    pdb_dir=indir+pdb+'/'
    pdb_dir_contents=os.listdir(pdb_dir)
    for item in pdb_dir_contents:
        if os.path.isdir(pdb_dir+item) and int(item)>num_models_offset and int(item)<=num_models_offset+num_models_per_PDB:
            model_subdir=pdb_dir+item+'/'
            model_subdir_contents=os.listdir(model_subdir)
            for item2 in model_subdir_contents:
                if '.o' in item2:
                    stats=[]
                    text=functions_lib.run_return('head -n 4 '+model_subdir+item2)
                    text+=functions_lib.run_return('tail -n 12 '+model_subdir+item2)
                    lines=text.split('\n')
                    #print text
                    if lines[-3]=='end_date:':
                        start_time=int(datetime.datetime.strftime(datetime.datetime.strptime(lines[3],input_time_format),output_time_format))
                        end_time=int(datetime.datetime.strftime(datetime.datetime.strptime(lines[-2],input_time_format),output_time_format))
                        total_energy=None
                        loop_rms=None
                        for line in lines:
                            if 'total_energy' in line:
                                total_energy=float(line.split('total_energy')[1].strip(': '))
                            elif 'loop_rms' in line:
                                loop_rms=float(line.split('loop_rms')[1].strip(': '))
                        #--
                        if total_energy!=None and loop_rms!=None:
                            #create new model
                            model=Model()
                            model.id=item
                            model.loop_rms=loop_rms
                            model.total_energy=total_energy
                            model.runtime=end_time-start_time
                            models.append(model)
                            runtimes.append(model.runtime)
                            #write to outfile
                            outstring=pdb+'\t'+model.id+'\t'+str(model.loop_rms)+'\t'+str(model.total_energy)+'\t'+str(model.runtime)
                            #print outstring
                            outfile.write(outstring+'\n')
    #------
    print len(models),'models'
    lowest_energy_model=Model()
    lowest_rmsd_model=Model()
    for model in models:
        if model.total_energy<lowest_energy_model.total_energy:
            lowest_energy_model=model
        #-
        if model.loop_rms<lowest_rmsd_model.loop_rms:
            lowest_rmsd_model=model
    #--
    print 'best model (i.e. lowest energy):',lowest_energy_model.id,lowest_energy_model.total_energy,lowest_energy_model.loop_rms
    print 'closest model (i.e. lowest rmsd):',lowest_rmsd_model.id,lowest_rmsd_model.total_energy,lowest_rmsd_model.loop_rms
    best_models.append(lowest_energy_model)
    closest_models.append(lowest_rmsd_model)
#-
print
print 'Global statistics ((min, median and max) energy and rmsd):'
print 'best models:',getEnergyStats(best_models),getRmsdStats(best_models)
print 'closest models:',getEnergyStats(closest_models),getRmsdStats(closest_models)
print 'average runtime:',Statistics.averageAndStddev(runtimes)
outfile.close()
print
print outfile_name                


end_time=time.time()
print "\ntime consumed: "+str(end_time-start_time)
