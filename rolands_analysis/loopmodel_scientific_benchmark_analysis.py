#!/usr/bin/python
#This file was developed and written by Roland A. Pache, Copyright (C) 2011, 2012.

#import libraries
import re
import os
import time
import sys
import math
from subprocess import *
#from optparse import OptionParser
#import Parameters
import Statistics
import functions_lib
#import Proteins


print '\nAnalyses the loopmodel scientific benchmark.\n'


#constants
dummy_eps='/home/rpache/postdoc/projects/KIC_with_fragments/results/dummy_image.eps'
top_X=5
#['refactored_KIC_max_2k','refactored_KIC_max_10k','KIC_with_fragments_nohoms_max_2k','KIC_with_fragments_nohoms_max_10k','KIC_with_fragments_nohoms']
#['CCD','CCD_nohoms','CCD_nohoms_ramaswitch','legacy_KIC','refactored_KIC','refactored_KIC_bugfix_linmin','refactored_KIC_bugfix_dfpmin','next_generation_KIC','KIC_with_fragments','KIC_with_fragments_nohoms']
#['CCD','CCD_nohoms','legacy_KIC','refactored_KIC','next_generation_KIC','KIC_with_fragments','KIC_with_fragments_nohoms','perturb_KIC_refactor_refine_KIC_with_fragments_nohoms','perturb_KIC_with_fragments_nohoms_refine_KIC_refactor']
#['CCD','CCD_nohoms','legacy_KIC','refactored_KIC','next_generation_KIC','KIC_with_fragments','KIC_with_fragments_nohoms']
#['CCD','CCD_nohoms','legacy_KIC','refactored_KIC','next_generation_KIC','KIC_with_fragments','KIC_with_fragments_nohoms','perturb_KIC_refactor_refine_KIC_with_fragments_nohoms','perturb_KIC_with_fragments_nohoms_refine_KIC_refactor']
#['CCD','CCD_nohoms','legacy_KIC','refactored_KIC','next_generation_KIC','KIC_with_fragments','KIC_with_fragments_nohoms']
#['CCD','CCD_nohoms','legacy_KIC','refactored_KIC','next_generation_KIC','next_generation_KIC_only_3_outer_cycles','KIC_with_fragments','KIC_with_fragments_nohoms','perturb_KIC_refactor_refine_KIC_with_fragments_nohoms','perturb_KIC_with_fragments_nohoms_refine_KIC_refactor']
#['next_generation_KIC','KIC_with_fragments_nohoms']
#['CCD','CCD_nohoms','legacy_KIC','refactored_KIC','next_generation_KIC','KIC_with_fragments','KIC_with_fragments_nohoms','perturb_KIC_refactor_refine_KIC_with_fragments_nohoms','perturb_KIC_with_fragments_nohoms_refine_KIC_refactor']

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
    return (min(energies),Statistics.median(energies),Statistics.average(energies),max(energies))

def getRmsdStats(models):
    rmsds=[]
    for model in models:
        if model.id!=None:
            rmsds.append(model.loop_rms)
    #--
    return (min(rmsds),Statistics.median(rmsds),Statistics.average(rmsds),max(rmsds))

def getRuntimeStats(models):
    runtimes=[]
    for model in models:
        if model.id!=None:
            runtimes.append(model.runtime)
    #--
    return (min(runtimes),Statistics.median(runtimes),Statistics.average(runtimes),max(runtimes))

def scatterplot(sorted_models,input_structure_score,rgb_color,outfile_name):
    outfile=open(outfile_name,'w')
    outfile.write('#Model\tLoop_rmsd\tTotal_energy\n')
    outfile.write('input_structure\t0.0\t'+str(input_structure_score)+'\n')
    outfile.write('\n\n')
    energies=[]
    num_subangstrom_models=0
    for model in sorted_models:
        outfile.write(model.id+'\t'+str(model.loop_rms)+'\t'+str(model.total_energy)+'\n')
        energies.append(model.total_energy)
        #check if model is sub-Angstrom
        if model.loop_rms<1.0:
            num_subangstrom_models+=1
        #-
    #-
    percent_subangstrom_models=round(100*num_subangstrom_models/float(len(sorted_models)),2)
    min_energy=min(energies)
    max_energy=max(energies)
    third_quartile=Statistics.quantile(energies,0.75)
    outfile.write('\n\n')
    i=4
    while i>=0:
        model=sorted_models[i]
        outfile.write(model.id+'\t'+str(model.loop_rms)+'\t'+str(model.total_energy)+'\n')
        i-=1
    #-
    outfile.write('\n\n')
    best_model=sorted_models[0]
    #when looking for the best model, consider the top X lowest energy models and pick the one with lowest rmsd
    for i in range(top_X):
        best_model_candidate=sorted_models[i]
        if best_model_candidate.loop_rms<best_model.loop_rms:
            best_model=best_model_candidate
    #--
    outfile.write(best_model.id+'\t'+str(best_model.loop_rms)+'\t'+str(best_model.total_energy)+'\n')
    outfile.close()
    pdb=outfile_name.split('/')[-1].split('_')[0]
    gnuplot_commands='\nset autoscale\
    \nset border 31\
    \nset tics out\
    \nset terminal postscript eps enhanced color "Helvetica" 24\
    \n#set size 1,1.5\
    \n#set size ratio 1\
    \n#set xtics ("default" 1, "default" 2, "H/Y" 3, "Y/H" 4, "default" 6, "default" 7, "H/Y" 8, "Y/H" 9) rotate by -45\
    \nset xtics autofreq\
    \nset xtics nomirror\
    \nset ytics autofreq\
    \nset ytics nomirror\
    \nset noy2tics\
    \nset nox2tics\
    \n\
    \nset style line 1 lt 1 lc rgb "dark-magenta" lw 2\
    \nset style line 2 lt 1 lc rgb "blue" lw 2 ps 1 pt 7\
    \nset style line 3 lt 1 lc rgb "forest-green" lw 2 ps 2 pt 13\
    \nset style line 4 lt 1 lc rgb "gold" lw 2 ps 1 pt 7\
    \nset style line 5 lt 1 lc rgb "red" lw 2 ps 2 pt 13\
    \nset style line 6 lt 1 lc rgb "black" lw 2\
    \nset style line 7 lt 1 lc rgb "dark-gray" lw 2\
    \nset style line 8 lt 1 lc rgb "gray" lw 2\
    \nset style line 9 lt 2 lc rgb "dark-gray" lw 5\
    \n\
    \nset boxwidth 0.75\
    \n\
    \nset key below right\
    \nset xrange [0:]\
    \nset encoding iso_8859_1\
    \nset title "'+pdb+': '+str(percent_subangstrom_models)+'% sub-{/E \305} models"\
    \nset xlabel "r.m.s. deviation to crystal loop [{/E \305}]"\
    \nset arrow from 1, graph 0 to 1, graph 1 ls 9 nohead\
    \nset ylabel "Rosetta all-atom score"\
    \nset output "'+outfile_name.split('.')[0]+'_all.eps"\
    \n#plot "'+outfile_name+'" index 0 using ($2):($3) with points ls 3 title "starting structure" axes x1y1,\
    \nplot "'+outfile_name+'" index 1 using ($2):($3) with points ls 2 title "all models" axes x1y1,\
    "'+outfile_name+'" index 2 using ($2):($3) with points ls 4 title "5 lowest energy models" axes x1y1,\
    "'+outfile_name+'" index 3 using ($2):($3) with points ls 5 title "top 5 best model" axes x1y1\
    \n\
    \nset yrange [:'+str(third_quartile)+']\
    \nset output "'+outfile_name.split('.')[0]+'_third_quartile.eps"\
    \nset xrange [0:]\
    \n#plot "'+outfile_name+'" index 0 using ($2):($3) with points ls 3 title "starting structure" axes x1y1,\
    \nplot "'+outfile_name+'" index 1 using ($2):($3) with points ls 2 title "75% lowest-scoring models" axes x1y1,\
    "'+outfile_name+'" index 2 using ($2):($3) with points ls 4 title "5 lowest energy models" axes x1y1,\
    "'+outfile_name+'" index 3 using ($2):($3) with points ls 5 title "top 5 best model" axes x1y1\
    \n'
    gnuplot_scriptname=outfile_name.split('.')[0]+'.gnu'
    functions_lib.newFile(gnuplot_commands,gnuplot_scriptname)
    functions_lib.run('gnuplot '+gnuplot_scriptname)
    return percent_subangstrom_models


def plotHistogram(sorted_models,percent_subangstrom_models,rgb_color,outfile_name):
    model_rmsds=[]
    for model in sorted_models:
        model_rmsds.append(model.loop_rms)
    #-
    tuples=Statistics.histogram(model_rmsds,100)
    outfile=open(outfile_name,'w')
    outfile.write('#All models\n')
    outfile.write('#RMSD\tFrequency\n')
    for (rmsd,count) in tuples:
        outfile.write(str(rmsd)+'\t'+str(100*count/float(len(model_rmsds)))+'\n')
    #-
    outfile.close()
    gnuplot_commands='\nset autoscale\
    \nset border 31\
    \nset tics out\
    \nset terminal postscript eps enhanced color "Helvetica" 24\
    \n#set size 1,1.5\
    \n#set size ratio 1\
    \n#set xtics ("default" 1, "default" 2, "H/Y" 3, "Y/H" 4, "default" 6, "default" 7, "H/Y" 8, "Y/H" 9) rotate by -45\
    \nset xtics autofreq\
    \nset xtics nomirror\
    \nset ytics autofreq\
    \nset ytics nomirror\
    \nset noy2tics\
    \nset nox2tics\
    \n\
    \nset style line 1 lt 1 lc rgb "dark-magenta" lw 2\
    \nset style line 2 lt 1 lc rgb "blue" lw 8 ps 1 pt 7\
    \nset style line 3 lt 1 lc rgb "forest-green" lw 2 ps 2 pt 13\
    \nset style line 4 lt 1 lc rgb "gold" lw 2 ps 1 pt 7\
    \nset style line 5 lt 1 lc rgb "red" lw 2 ps 2 pt 13\
    \nset style line 6 lt 1 lc rgb "black" lw 2\
    \nset style line 7 lt 1 lc rgb "dark-gray" lw 2\
    \nset style line 8 lt 1 lc rgb "gray" lw 2\
    \nset style line 9 lt 2 lc rgb "dark-gray" lw 5\
    \n\
    \nset boxwidth 0.75\
    \n\
    \nset key below right\
    \nset xrange [0:]\
    \nset encoding iso_8859_1\
    \nset title "'+pdb+': '+str(percent_subangstrom_models)+'% sub-{/E \305} models"\
    \nset xlabel "r.m.s. deviation to crystal loop [{/E \305}]"\
    \nset yrange [0:]\
    \nset arrow from 1, graph 0 to 1, graph 1 ls 9 nohead\
    \nset ylabel "Fraction of models [%]"\
    \nset output "'+outfile_name.split('.')[0]+'.eps"\
    \nplot "'+outfile_name+'" index 0 using ($1):($2) smooth bezier with lines ls 2 title "all models" axes x1y1\
    \n'
    gnuplot_scriptname=outfile_name.split('.')[0]+'.gnu'
    functions_lib.newFile(gnuplot_commands,gnuplot_scriptname)
    functions_lib.run('gnuplot '+gnuplot_scriptname)
    return model_rmsds


def boxplots(models,percent_subangstrom_models_list,outfile_name):
    boxplot_data={}
    boxplot_data[1]=[]
    boxplot_data[2]=[]
    for model in models:
        boxplot_data[1].append(model.loop_rms)
        boxplot_data[2].append(model.total_energy)
    #-
    boxplot_data[3]=percent_subangstrom_models_list
    boxplot_results=Statistics.tukeyBoxAndWhisker(boxplot_data)
    outfile=open(outfile_name,'w')
    for x in boxplot_results:
        outfile.write('#x\t'+'lower\t'+'first_quartile\t'+'median\t'+'third_quartile\t'+'upper\n')
        (tuple,outliers)=boxplot_results[x]
        for item in tuple:
            outfile.write(str(item)+'\t')
        #-
        outfile.write('\n\n\n')
        outfile.write('#x\t'+'outlier\n')
        if len(outliers)==0:
            outfile.write(str(x)+'\t?\n')
        else:
            for outlier in outliers:
                outfile.write(str(x)+'\t'+str(outlier)+'\n')
        #--
        outfile.write('\n\n')
    #-
    outfile.close()
    gnuplot_commands='\nset autoscale\
    \nset border 31\
    \nset tics out\
    \nset terminal postscript eps enhanced color "Helvetica" 24\
    \n#set size 1,1.5\
    \nset size ratio 1\
    \n#set xtics ("default" 1, "-correct" 2, "default" 3, "-correct" 4) rotate by -45\
    \nset noxtics\
    \nset xrange [0.5:1.5]\
    \nset nox2tics\
    \nset ytics 1\
    \nset ytics nomirror\
    \nset noy2tics\
    \n\
    \nset style line 1 lt 1 lc rgb "dark-magenta" lw 2\
    \nset style line 2 lt 1 lc rgb "blue" lw 5 pt 7\
    \nset style line 3 lt 1 lc rgb "forest-green" lw 5\
    \nset style line 4 lt 1 lc rgb "gold" lw 2\
    \nset style line 5 lt 1 lc rgb "red" lw 5 pt 7\
    \nset style line 6 lt 1 lc rgb "black" lw 5\
    \nset style line 7 lt 1 lc rgb "dark-gray" lw 2\
    \nset style line 8 lt 1 lc rgb "gray" lw 2\
    \nset style line 9 lt 0 lc rgb "black" lw 5\
    \n\
    \nset boxwidth 0.25\
    \n\
    \nset key tmargin\
    \nset title "Best models performance distribution"\
    \nset noxlabel\
    \nset style fill solid 0.5\
    \nset encoding iso_8859_1\
    \nset ylabel "r.m.s. deviation to crystal loop [{/E \305}]"\
    \nset output "'+outfile_name.split('.')[0]+'_rmsd.eps"\
    \nf(x)=1\
    \nplot "'+outfile_name+'" index 0 using 1:3:2:6:5 with candlesticks whiskerbars ls 2 notitle axes x1y1,\
    "'+outfile_name+'" index 0 using 1:4:4:4:4 with candlesticks ls 6 notitle,\
    "'+outfile_name+'" index 1 using 1:2 with points ls 2 pt 6 notitle,\
    f(x) with lines ls 9 notitle\
    \n\
    \nset ylabel "Rosetta all-atom score"\
    \nset xrange [1.5:2.5]\
    \nset ytics autofreq\
    \nset output "'+outfile_name.split('.')[0]+'_energy.eps"\
    \nplot "'+outfile_name+'" index 2 using 1:3:2:6:5 with candlesticks whiskerbars ls 5 notitle axes x1y1,\
    "'+outfile_name+'" index 2 using 1:4:4:4:4 with candlesticks ls 6 notitle,\
    "'+outfile_name+'" index 3 using 1:2 with points ls 5 pt 6 notitle\
    \n\
    \nset title "Protocol performance distribution"\
    \nset ylabel "Fraction sub-{/E \305} models [%]"\
    \nset xrange [2.5:3.5]\
    \nset ytics 10\
    \nset output "'+outfile_name.split('.')[0]+'_percent_sub-A_models.eps"\
    \nplot "'+outfile_name+'" index 4 using 1:3:2:6:5 with candlesticks whiskerbars ls 3 notitle axes x1y1,\
    "'+outfile_name+'" index 4 using 1:4:4:4:4 with candlesticks ls 6 notitle,\
    "'+outfile_name+'" index 5 using 1:2 with points ls 3 pt 6 notitle\
    \n'
    gnuplot_scriptname=outfile_name.split('.')[0]+'.gnu'
    functions_lib.newFile(gnuplot_commands,gnuplot_scriptname)
    functions_lib.run('gnuplot '+gnuplot_scriptname)


            
start_time=time.time()

#parse input parameters
if len(sys.argv)!=2:
    print
    print 'Usage: ./loopmodel_scientific_benchmark_analysis.py PARAMETER_FILE'
    print
    sys.exit()
#-
parameter_file=sys.argv[1]
parameters=functions_lib.parseParameterFile(parameter_file)
main_dir=parameters['global_main_dir']
num_models_offset=int(parameters['loopmodel_num_models_offset'])
num_models_per_PDB=int(parameters['loopmodel_num_models_per_PDB'])
run_identifiers=parameters['loopmodel_run_identifiers'].split(',')
input_structures_dir=main_dir+parameters['loopmodel_structures_indir']
indir_prefix=main_dir+parameters['loopmodel_outdir']
model_start_index=num_models_offset+1
model_end_index=num_models_offset+num_models_per_PDB


#parse Rosetta score of input structures
print
print 'parsing total scores of input structures...'
input_structure_scores={}
pdb_ids=[]
input_structures_dir_contents=os.listdir(input_structures_dir)
for item in input_structures_dir_contents:
    if item.endswith('.pdb.gz'):
        score_line=functions_lib.run_return('gunzip -c '+input_structures_dir+item+'|grep pose')
        score=float(score_line.split()[-1])
        pdb_id=item.split('_')[0]
        print pdb_id,score
        pdb_ids.append(pdb_id)
        input_structure_scores[pdb_id]=score
#--
sorted_pdb_ids=sorted(pdb_ids)


#process all given run identifiers
percentage_subangstrom_models_map={}
best_models_rmsds_map={}
lowest_energy_models_rmsds_map={}
sampling_density_distributions_map={}
runtimes_map={}
rgb_colors=functions_lib.gnuplotColorWheel(len(run_identifiers), saturation_adjustment = 2)
for i in range(len(run_identifiers)):
    run_identifier=run_identifiers[i]
    rgb_color=rgb_colors[i]
    infile_name=indir_prefix+run_identifier+'.results'
    if not os.path.isfile(infile_name):
        continue
    #-
    outdir_prefix=indir_prefix+run_identifier+'_analysis/'
    print
    print infile_name
    print outdir_prefix
    
    #parse models
    runtimes=[]
    models_per_pdb={}
    all_models=[]
    infile=open(infile_name)
    for line in infile:
        if not line.startswith('#'):
            data=line.strip('\n').split('\t')
            if len(data)>4:
                pdb=data[0]
                if pdb not in models_per_pdb:
                     models_per_pdb[pdb]=[]
                #-
                model_index=int(data[1])
                if model_index>=model_start_index and model_index<=model_end_index:
                    model=Model()
                    model.id=pdb+'_'+str(model_index)
                    model.loop_rms=float(data[2])
                    model.total_energy=float(data[3])
                    model.runtime=int(data[4])
                    models_per_pdb[pdb].append(model)
                    all_models.append(model)
                    runtimes.append(model.runtime)
    #----
    infile.close()
    runtimes_map[run_identifier]=runtimes
    print 'total runtime [hours]:',int(sum(runtimes)/float(3600))
    (average_runtime,runtime_stddev)=Statistics.averageAndStddev(runtimes)
    print 'average runtime [seconds]:',int(round(average_runtime,0)),'+/-',int(round(runtime_stddev,0))

    #compute basic statistics and create rmsd vs. Rosetta score plots per PDB
    tex_tables=[]
    best_models=[]
    lowest_energy_models=[]
    closest_models=[]
    percent_subangstrom_models_list=[]
    print len(sorted_pdb_ids),'PDBs'
    tex_table_string='\\begin{tabular}{rr|rrr|rrr}\n\
    PDB &\# models &Top '+str(top_X)+' best model &Loop rmsd &Energy &Closest model &Loop rmsd &Energy\\\\\\hline\n'
    for pdb in sorted_pdb_ids:
        print
        print pdb
        if pdb not in models_per_pdb:
            print 'WARNING: no models found for',pdb
            continue
        #-
        if pdb not in sampling_density_distributions_map:
            sampling_density_distributions_map[pdb]={}
        #-
        models=models_per_pdb[pdb]
        #create outdir
        outdir=outdir_prefix+pdb+'/'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #-
        print len(models),'successful models'
        #determine best and closest model for the given pdb
        sorted_energy_models=sorted(models,lambda x, y: cmp(x.total_energy,y.total_energy))
        sorted_rmsd_models=sorted(models,lambda x, y: cmp(x.loop_rms,y.loop_rms))
        lowest_energy_model=sorted_energy_models[0]
        lowest_energy_models.append(lowest_energy_model)
        lowest_rmsd_model=sorted_rmsd_models[0]
        #when looking for the best model, consider the top X lowest energy models and pick the one with lowest rmsd
        best_model=lowest_energy_model
        for i in range(top_X):
            best_model_candidate=sorted_energy_models[i]
            if best_model_candidate.loop_rms<best_model.loop_rms:
                best_model=best_model_candidate
        #--
        print 'lowest energy model:',lowest_energy_model.id,lowest_energy_model.loop_rms,lowest_energy_model.total_energy
        print 'best model (i.e. lowest rmsd of top '+str(top_X)+' lowest energy models):',best_model.id,best_model.loop_rms,best_model.total_energy
        print 'closest model (i.e. lowest rmsd):',lowest_rmsd_model.id,lowest_rmsd_model.loop_rms,lowest_rmsd_model.total_energy
        best_models.append(best_model)
        closest_models.append(lowest_rmsd_model)
        #create scatterplot for each pdb
        outfile_name=outdir+pdb+'_models.out'
        input_structure_score=input_structure_scores[pdb]
        percent_subangstrom_models=scatterplot(sorted_energy_models,input_structure_score,rgb_color,outfile_name)
        percent_subangstrom_models_list.append(percent_subangstrom_models)
        #create sampling histogram for each pdb
        outfile_name=outdir+pdb+'_sampling_histogram.out'
        model_rmsds=plotHistogram(sorted_rmsd_models,percent_subangstrom_models,rgb_color,outfile_name)
        sampling_density_distributions_map[pdb][run_identifier]=model_rmsds
        #write tex table results row
        tex_table_string+=pdb+' &'+str(len(models))+' &'+best_model.id.split('_')[-1]+' &'+'%0.2f' % best_model.loop_rms+' &'+'%0.2f' % best_model.total_energy+' &'+lowest_rmsd_model.id.split('_')[-1]+' &'+'%0.2f' % lowest_rmsd_model.loop_rms+' &'+'%0.2f' % lowest_rmsd_model.total_energy+'\\\\\n'
    #-
    tex_table_string+='\\end{tabular}\n'
    tex_outfile_name=outdir_prefix+'individual_results.tex'
    tex_tables.append(tex_outfile_name)
    tex_outfile=open(tex_outfile_name,'w')
    tex_outfile.write(tex_table_string)
    tex_outfile.close()
    overall_stats=(best_models,closest_models)
    best_models_rmsds_map[run_identifier]=best_models
    lowest_energy_models_rmsds_map[run_identifier]=lowest_energy_models
    print
    print tex_outfile_name
    
    #create rmsd and energy boxplots for the best models
    outdir=outdir_prefix
    outfile_name=outdir+'best_models_dists.out'
    boxplots(best_models,percent_subangstrom_models_list,outfile_name)
    print outfile_name

    #calculate global stats across all pdbs and write overall performance tex table
    print
    print 'Global statistics (median rmsd, energy and runtime):'
    outfile_name=outdir_prefix+'overall_results.tex'
    tex_tables.append(outfile_name)
    percentage_subangstrom_models_map[run_identifier]=percent_subangstrom_models_list
    outfile=open(outfile_name,'w')
    outfile.write('\
    \n{\\bf Median fraction of sub-\\AA~models: '+str(round(Statistics.median(percent_subangstrom_models_list),2))+'\\%}\\\\[0.5cm]\
    \n\
    \\begin{tabular}{lrrr}\
    \nModel selection &Median loop rmsd &Median energy &Median runtime\\\\\\hline\
    \n')
    best_models_median_energy=round(getEnergyStats(best_models)[1],2)
    best_models_median_rmsd=round(getRmsdStats(best_models)[1],2)
    best_models_median_runtime=round(getRuntimeStats(best_models)[1],2)
    lowest_energy_models_median_energy=round(getEnergyStats(lowest_energy_models)[1],2)
    lowest_energy_models_median_rmsd=round(getRmsdStats(lowest_energy_models)[1],2)
    lowest_energy_models_median_runtime=round(getRuntimeStats(lowest_energy_models)[1],2)
    closest_models_median_energy=round(getEnergyStats(closest_models)[1],2)
    closest_models_median_rmsd=round(getRmsdStats(closest_models)[1],2)
    closest_models_median_runtime=round(getRuntimeStats(closest_models)[1],2)
    all_models_median_energy=round(getEnergyStats(all_models)[1],2)
    all_models_median_rmsd=round(getRmsdStats(all_models)[1],2)
    all_models_median_runtime=round(getRuntimeStats(all_models)[1],2)
    print 'best models median rmsd, energy and runtime:',best_models_median_rmsd,best_models_median_energy,best_models_median_runtime
    print 'lowest energy models median rmsd, energy and runtime:',lowest_energy_models_median_energy,lowest_energy_models_median_rmsd,lowest_energy_models_median_runtime
    print 'closest models median rmsd, energy and runtime:',closest_models_median_rmsd,closest_models_median_energy,closest_models_median_runtime
    print 'all models median rmsd, energy and runtime:',all_models_median_rmsd,all_models_median_energy,all_models_median_runtime
    outfile.write('{\\bf Top '+str(top_X)+' best model} &{\\bf '+'%0.2f' % best_models_median_rmsd+'} &{\\bf '+'%0.2f' % best_models_median_energy+'} &{\\bf '+'%0.2f' % best_models_median_runtime+'}\\\\\n')
    outfile.write('{\\bf Lowest energy model} &{\\bf '+'%0.2f' % lowest_energy_models_median_rmsd+'} &{\\bf '+'%0.2f' % lowest_energy_models_median_energy+'} &{\\bf '+'%0.2f' % lowest_energy_models_median_runtime+'}\\\\\n')
    outfile.write('Closest model &'+'%0.2f' % closest_models_median_rmsd+' &'+'%0.2f' % closest_models_median_energy+' &'+'%0.2f' % closest_models_median_runtime+'\\\\\n')
    outfile.write('{\\bf All models} &{\\bf '+'%0.2f' % all_models_median_rmsd+'} &{\\bf '+'%0.2f' % all_models_median_energy+'} &{\\bf '+'%0.2f' % all_models_median_runtime+'}\\\\\n')
    outfile.write('\\end{tabular}\n')
    outfile.close()
    print outfile_name

    #put all model output figures into a tex table
    num_pdbs=len(sorted_pdb_ids)
    num_rows=5
    num_cols=2
    num_plots_per_page=num_rows*num_cols
    num_pages=int(math.ceil(2*num_pdbs/float(num_plots_per_page)))
    print num_pdbs,'pdbs'
    print num_pages,'pages'
    index=0
    for i in range(num_pages):
        outfile1_name=outdir_prefix+'/all_models_'+str(i+1)+'.tex'
        outfile2_name=outdir_prefix+'/third_quartile_models_'+str(i+1)+'.tex'
        outstring1='\\begin{figure}\n\
        \\vspace*{-2cm}\n\
        \\resizebox{\\textwidth}{!}{\n\
        \\begin{tabular}{'
        for j in range(num_cols):
            outstring1+='c'
        #-
        outstring1+='}\n'
        outstring2=outstring1
        for j in range(num_rows):
            k=0
            while k<num_cols:
                if index<num_pdbs:
                    pdb=sorted_pdb_ids[index]
                    if pdb in models_per_pdb:
                        outstring1+='\includegraphics{'+outdir_prefix+pdb+'/'+pdb+'_models_all.eps} &'+'\includegraphics{'+outdir_prefix+pdb+'/'+pdb+'_sampling_histogram.eps} &'
                        outstring2+='\includegraphics{'+outdir_prefix+pdb+'/'+pdb+'_models_third_quartile.eps} &'+'\includegraphics{'+outdir_prefix+pdb+'/'+pdb+'_sampling_histogram.eps} &'
                    else:
                        outstring1+='\includegraphics{'+dummy_eps+'} &'+'\includegraphics{'+dummy_eps+'} &'
                        outstring2+='\includegraphics{'+dummy_eps+'} &'+'\includegraphics{'+dummy_eps+'} &'
                #--
                else:
                    outstring1+='\includegraphics{'+dummy_eps+'} &'+'\includegraphics{'+dummy_eps+'} &'
                    outstring2+='\includegraphics{'+dummy_eps+'} &'+'\includegraphics{'+dummy_eps+'} &'
                #-
                k+=2
                index+=1
            #-
            outstring1=outstring1.rstrip(' &')+'\\\\\n'
            outstring2=outstring2.rstrip(' &')+'\\\\\n'
        #-
        outstring1+='\\end{tabular}\n\
        }\\end{figure}\n'
        outstring2+='\\end{tabular}\n\
        }\\end{figure}\n'
        outfile=open(outfile1_name,'w')
        outfile.write(outstring1)
        outfile.close()
        outfile=open(outfile2_name,'w')
        outfile.write(outstring2)
        outfile.close()
        print outfile1_name
        print outfile2_name
    #-
#-


#create protocols comparison bar chart and box plots
print
print 'creating protocols comparison figures...'
outdir_prefix=indir_prefix+'protocols_comparison/'
if not os.path.isdir(outdir_prefix):
    os.makedirs(outdir_prefix)
#-
outfile1_name=outdir_prefix+'percentage_subangstrom_models.out'
outfile=open(outfile1_name,'w')
reversed_run_identifiers=list(reversed(run_identifiers))
reversed_rgb_colors=list(reversed(rgb_colors))
for i in range(len(reversed_run_identifiers)):
    outfile.write('#Protocol\tx\tlower\tfirst_quartile\tmedian\tthird_quartile\tupper\n')
    run_identifier=reversed_run_identifiers[i]
    ((x,lower,first_quartile,median,third_quartile,upper),outliers)=(('?','?','?','?','?','?'),[])
    if run_identifier in percentage_subangstrom_models_map:
        boxplot_data={}
        boxplot_data[1]=percentage_subangstrom_models_map[run_identifier]
        boxplot_results=Statistics.tukeyBoxAndWhisker(boxplot_data)
        ((x,lower,first_quartile,median,third_quartile,upper),outliers)=boxplot_results[1]
    #-
    outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(lower)+'\t'+str(first_quartile)+'\t'+str(median)+'\t'+str(third_quartile)+'\t'+str(upper)+'\n\n\n')
    outfile.write('#Protocol\tx\toutlier\n')
    if len(outliers)==0:
        outfile.write(run_identifier+'\t'+str(i+1)+'\t?\n')
    else:
        for outlier in outliers:
            outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(outlier)+'\n')
    #--
    outfile.write('\n\n')
#-
outfile.close()
outfile2_name=outdir_prefix+'best_models_rmsd.out'
outfile=open(outfile2_name,'w')
for i in range(len(reversed_run_identifiers)):
    outfile.write('#Protocol\tx\tlower\tfirst_quartile\tmedian\tthird_quartile\tupper\n')
    run_identifier=reversed_run_identifiers[i]
    ((x,lower,first_quartile,median,third_quartile,upper),outliers)=(('?','?','?','?','?','?'),[])
    if run_identifier in best_models_rmsds_map:
        best_models=best_models_rmsds_map[run_identifier]
        boxplot_data={}
        boxplot_data[1]=[]
        for model in best_models:
            boxplot_data[1].append(model.loop_rms)
        #-
        boxplot_results=Statistics.tukeyBoxAndWhisker(boxplot_data)
        ((x,lower,first_quartile,median,third_quartile,upper),outliers)=boxplot_results[1]
    #-
    outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(lower)+'\t'+str(first_quartile)+'\t'+str(median)+'\t'+str(third_quartile)+'\t'+str(upper)+'\n\n\n')
    outfile.write('#Protocol\tx\toutlier\n')
    if len(outliers)==0:
        outfile.write(run_identifier+'\t'+str(i+1)+'\t?\n')
    else:
        for outlier in outliers:
            outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(outlier)+'\n')
    #--
    outfile.write('\n\n')        
#-
outfile.close()
outfile3_name=outdir_prefix+'lowest_energy_models_rmsd.out'
outfile=open(outfile3_name,'w')
for i in range(len(reversed_run_identifiers)):
    outfile.write('#Protocol\tx\tlower\tfirst_quartile\tmedian\tthird_quartile\tupper\n')
    run_identifier=reversed_run_identifiers[i]
    ((x,lower,first_quartile,median,third_quartile,upper),outliers)=(('?','?','?','?','?','?'),[])
    if run_identifier in lowest_energy_models_rmsds_map:
        lowest_energy_models=lowest_energy_models_rmsds_map[run_identifier]
        boxplot_data={}
        boxplot_data[1]=[]
        for model in lowest_energy_models:
            boxplot_data[1].append(model.loop_rms)
        #-
        boxplot_results=Statistics.tukeyBoxAndWhisker(boxplot_data)
        ((x,lower,first_quartile,median,third_quartile,upper),outliers)=boxplot_results[1]
    #-
    outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(lower)+'\t'+str(first_quartile)+'\t'+str(median)+'\t'+str(third_quartile)+'\t'+str(upper)+'\n\n\n')
    outfile.write('#Protocol\tx\toutlier\n')
    if len(outliers)==0:
        outfile.write(run_identifier+'\t'+str(i+1)+'\t?\n')
    else:
        for outlier in outliers:
            outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(outlier)+'\n')
    #--
    outfile.write('\n\n')
#-
outfile.close()
outfile4_name=outdir_prefix+'all_models_runtime.out'
outfile=open(outfile4_name,'w')
for i in range(len(reversed_run_identifiers)):
    outfile.write('#Protocol\tx\tlower\tfirst_quartile\tmedian\tthird_quartile\tupper\n')
    run_identifier=reversed_run_identifiers[i]
    ((x,lower,first_quartile,median,third_quartile,upper),outliers)=(('?','?','?','?','?','?'),[])
    if run_identifier in lowest_energy_models_rmsds_map:
        boxplot_data={}
        boxplot_data[1]=[]
        runtimes=runtimes_map[run_identifier]
        for runtime in runtimes:
            boxplot_data[1].append(runtime/float(60))
        #-
        boxplot_results=Statistics.tukeyBoxAndWhisker(boxplot_data)
        ((x,lower,first_quartile,median,third_quartile,upper),outliers)=boxplot_results[1]
    #-
    outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(lower)+'\t'+str(first_quartile)+'\t'+str(median)+'\t'+str(third_quartile)+'\t'+str(upper)+'\n\n\n')
    outfile.write('#Protocol\tx\toutlier\n')
    if len(outliers)==0:
        outfile.write(run_identifier+'\t'+str(i+1)+'\t?\n')
    else:
        for outlier in outliers:
            outfile.write(run_identifier+'\t'+str(i+1)+'\t'+str(outlier)+'\n')
    #--
    outfile.write('\n\n')
#-
outfile.close()
gnuplot_commands='\nset autoscale\
\nset border 31\
\nset tics out\
\nset terminal postscript eps enhanced color "Helvetica" 24\
\nset size 1,3\
\n#set size ratio 0.5\
\nset xtics ("'
for i in range(len(reversed_run_identifiers)):
    run_identifier=reversed_run_identifiers[i]
    if '2k' in run_identifier or '10k' in run_identifier:
        gnuplot_commands+=run_identifier.replace('_',' ')+'" '+str(i+1)+', "'
    else:
        gnuplot_commands+=run_identifier.replace('_',' ')+' max 100k" '+str(i+1)+', "'
#-
gnuplot_commands=gnuplot_commands.rstrip(', "')+') rotate by -90\
\nset xtics nomirror\
\nset ytics autofreq rotate by -90\
\nset ytics nomirror\
\nset noy2tics\
\nset nox2tics\
\n\
\nset style line 1 lt 1 lc rgb "dark-magenta" lw 2\
\nset style line 2 lt 1 lc rgb "blue" lw 5 ps 1 pt 7\
\nset style line 3 lt 1 lc rgb "forest-green" lw 2 ps 2 pt 13\
\nset style line 4 lt 1 lc rgb "gold" lw 2 ps 1 pt 7\
\nset style line 5 lt 1 lc rgb "red" lw 2 ps 2 pt 13\
\nset style line 6 lt 1 lc rgb "black" lw 2\
\nset style line 7 lt 1 lc rgb "dark-gray" lw 2\
\nset style line 8 lt 1 lc rgb "gray" lw 2\
\nset style line 9 lt 2 lc rgb "dark-gray" lw 5\
\n\
\nset boxwidth 0.75\
\n\
\nset key below right\
\nset xrange [0:'+str(len(reversed_run_identifiers)+1)+']\
\nset encoding iso_8859_1\
\n#set title "Sub-{/E \305} sampling performance (based on all models)"\
\nset notitle\
\nunset xlabel\
\nset yrange [0:]\
\nset ylabel "Median fraction of sub-{/E \305} models [%]" rotate by -90\
\nset output "'+outfile1_name.split('.')[0]+'_barchart.eps"\
\nplot '
for i in range(len(reversed_run_identifiers)):
    gnuplot_commands+='"'+outfile1_name+'" index '+str(2*i)+' using 2:5 with boxes fs solid 0.5 lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 notitle, '
#-
gnuplot_commands=gnuplot_commands.rstrip(', ')+'\n\
\nset ylabel "Fraction of sub-{/E \305} models [%]" rotate by -90\
\nset output "'+outfile1_name.split('.')[0]+'.eps"\
\nset style fill solid 0.5\
\nset ytics autofreq rotate by -45\
\nplot '
for i in range(len(reversed_run_identifiers)):
    gnuplot_commands+='"'+outfile1_name+'" index '+str(2*i)+' using 2:4:3:7:6 with candlesticks whiskerbars lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 notitle,\
    "'+outfile1_name+'" index '+str(2*i)+' using 2:5:5:5:5 with candlesticks lt 1 lc rgb "black" lw 5 notitle,\
    "'+outfile1_name+'" index '+str(2*i+1)+' using 2:3 with points lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 pt 6 notitle,'
#-
gnuplot_commands=gnuplot_commands.rstrip(', ')+'\n\
\n#set title "Loop reconstruction performance (based on top '+str(top_X)+' models)"\
\nset ylabel "r.m.s. deviation to crystal loop [{/E \305}]" rotate by -90\
\nset ytics autofreq rotate by -90\
\nset output "'+outfile2_name.split('.')[0]+'.eps"\
\nf(x)=1\
\nplot '
for i in range(len(reversed_run_identifiers)):
    gnuplot_commands+='"'+outfile2_name+'" index '+str(2*i)+' using 2:4:3:7:6 with candlesticks whiskerbars lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 notitle,\
    "'+outfile2_name+'" index '+str(2*i)+' using 2:5:5:5:5 with candlesticks lt 1 lc rgb "black" lw 5 notitle,\
    "'+outfile2_name+'" index '+str(2*i+1)+' using 2:3 with points lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 pt 6 notitle,'
#-
gnuplot_commands+='\
f(x) with lines ls 9 notitle\n\
\n#set title "Loop reconstruction performance (based on lowest energy models)"\
\nset output "'+outfile3_name.split('.')[0]+'.eps"\
\nplot '
for i in range(len(reversed_run_identifiers)):
    gnuplot_commands+='"'+outfile3_name+'" index '+str(2*i)+' using 2:4:3:7:6 with candlesticks whiskerbars lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 notitle,\
    "'+outfile3_name+'" index '+str(2*i)+' using 2:5:5:5:5 with candlesticks lt 1 lc rgb "black" lw 5 notitle,\
    "'+outfile3_name+'" index '+str(2*i+1)+' using 2:3 with points lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 pt 6 notitle,'
#-
gnuplot_commands+='\
f(x) with lines ls 9 notitle\n\
\n#set title "Runtime (based on all models)"\
\nset ylabel "Runtime [min]" rotate by -90\
\nset ytics 30 rotate by -45\
\nset output "'+outfile4_name.split('.')[0]+'.eps"\
\nplot '
for i in range(len(reversed_run_identifiers)):
    gnuplot_commands+='"'+outfile4_name+'" index '+str(2*i)+' using 2:4:3:7:6 with candlesticks whiskerbars lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 notitle,\
    "'+outfile4_name+'" index '+str(2*i)+' using 2:5:5:5:5 with candlesticks lt 1 lc rgb "black" lw 5 notitle,\
    "'+outfile4_name+'" index '+str(2*i+1)+' using 2:3 with points lt 1 lc rgb "#'+reversed_rgb_colors[i]+'" lw 5 pt 6 notitle,'
#-
gnuplot_scriptname=outdir_prefix+'performance_comparison.gnu'
functions_lib.newFile(gnuplot_commands.rstrip(', '),gnuplot_scriptname)
functions_lib.run('gnuplot '+gnuplot_scriptname)


#create protocols comparison density distributions
for pdb in sorted_pdb_ids:
    outdir=outdir_prefix+pdb+'/'
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #-
    outfile_name=outdir+pdb+'_sampling_density_distributions.out'
    outfile=open(outfile_name,'w')
    outfile.write('#Protocol\tRMSD\tFrequency\n')
    for i in range(len(run_identifiers)):
        run_identifier=run_identifiers[i]
        if run_identifier not in sampling_density_distributions_map[pdb]:
            outfile.write(run_identifier+'\t?\t?\n')
        else:
            model_rmsds=sampling_density_distributions_map[pdb][run_identifier]
            tuples=Statistics.histogram(model_rmsds,100)
            for (rmsd,count) in tuples:
                outfile.write(run_identifier+'\t'+str(rmsd)+'\t'+str(100*count/float(len(model_rmsds)))+'\n')
        #--
        outfile.write('\n\n')
    #-
    outfile.close()
    gnuplot_commands='\nset autoscale\
    \nset border 31\
    \nset tics out\
    \nset terminal postscript eps enhanced color "Helvetica" 24\
    \n#set size 1.5,1\
    \n#set size ratio 1\
    \n#set xtics ("default" 1, "default" 2, "H/Y" 3, "Y/H" 4, "default" 6, "default" 7, "H/Y" 8, "Y/H" 9) rotate by -45\
    \nset xtics autofreq\
    \nset xtics nomirror\
    \nset ytics autofreq\
    \nset ytics nomirror\
    \nset noy2tics\
    \nset nox2tics\
    \n\
    \nset style line 1 lt 1 lc rgb "dark-magenta" lw 2\
    \nset style line 2 lt 1 lc rgb "blue" lw 5 ps 1 pt 7\
    \nset style line 3 lt 1 lc rgb "forest-green" lw 2 ps 2 pt 13\
    \nset style line 4 lt 1 lc rgb "gold" lw 2 ps 1 pt 7\
    \nset style line 5 lt 1 lc rgb "red" lw 2 ps 2 pt 13\
    \nset style line 6 lt 1 lc rgb "black" lw 2\
    \nset style line 7 lt 1 lc rgb "dark-gray" lw 2\
    \nset style line 8 lt 1 lc rgb "gray" lw 2\
    \nset style line 9 lt 2 lc rgb "dark-gray" lw 5\
    \n\
    \nset boxwidth 0.75\
    \n\
    \n#set key below right\
    \nunset key\
    \nset xrange [0:]\
    \nset encoding iso_8859_1\
    \nset title "'+pdb+'"\
    \nset xlabel "r.m.s. deviation to crystal loop [{/E \305}]"\
    \nset yrange [0:]\
    \nset arrow from 1, graph 0 to 1, graph 1 ls 9 nohead\
    \nset ylabel "Fraction of models [%]"\
    \nset output "'+outfile_name.split('.')[0]+'.eps"\
    \nplot '
    for i in range(len(run_identifiers)):
        run_identifier=run_identifiers[i]
        gnuplot_commands+='"'+outfile_name+'" index '+str(i)+' using 2:3 smooth bezier with lines lt 1 lc rgb "#'+rgb_colors[i]+'" lw 8 title "'+run_identifier.replace('_',' ')+'", '
    #-
    gnuplot_scriptname=outfile_name.split('.')[0]+'.gnu'
    functions_lib.newFile(gnuplot_commands.rstrip(', '),gnuplot_scriptname)
    functions_lib.run('gnuplot '+gnuplot_scriptname)
#-
#put all model output figures into a tex table
num_pdbs=len(sorted_pdb_ids)
num_rows=3
num_cols=3
num_plots_per_page=num_rows*num_cols
num_pages=int(math.ceil(num_pdbs/float(num_plots_per_page)))
print num_pdbs,'pdbs'
print num_pages,'pages'
index=0
for i in range(num_pages):
    outfile_name=outdir_prefix+'sampling_density_distributions_'+str(i+1)+'.tex'
    outstring='\\begin{sidewaysfigure}\n\
    \\resizebox{\\textwidth}{!}{\n\
    \\begin{tabular}{'
    for j in range(num_cols):
        outstring+='c'
    #-
    outstring+='}\n'
    for j in range(num_rows):
        k=0
        while k<num_cols:
            if index<num_pdbs:
                pdb=sorted_pdb_ids[index]
                if pdb in sampling_density_distributions_map:
                    outstring+='\includegraphics{'+outdir_prefix+pdb+'/'+pdb+'_sampling_density_distributions.eps} &'
                else:
                    outstring+='\includegraphics{'+dummy_eps+'} &'
            #--
            else:
                outstring+='\includegraphics{'+dummy_eps+'} &'
            #-
            k+=1
            index+=1
        #-
        outstring=outstring.rstrip(' &')+'\\\\\n'
    #-
    outstring+='\\end{tabular}\n\
    }\\end{sidewaysfigure}\n'
    outfile=open(outfile_name,'w')
    outfile.write(outstring)
    outfile.close()
    print outfile_name
#-
print
print outdir_prefix
    

end_time=time.time()
print
print "\ntime consumed: "+str(end_time-start_time)
