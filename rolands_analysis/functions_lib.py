#This file was developed and written by Roland A. Pache, Ph.D., Copyright (C) 2011, 2012.

import os
import shutil
import math
import colorsys
from subprocess import *


#creates the directory if it does not yet exist, and deletes its contents otherwise
#input: name of the directory
#output: -
def createDir(dir):
    if os.path.isdir(dir):
        shutil.rmtree(dir)
    #-
    os.makedirs(dir)


#parser for parameter file
#input: name of the parameter file
#output: map: parameter -> value
def parseParameterFile(param_filename):
    parameters={}
    file=open(param_filename)
    lines=file.read().split('\n')
    file.close()
    for line in lines:
        data=line.split('=')
        if len(data)>1:
            #parameter
            parameter=data[0].lstrip(' ').rstrip(' ')
            #value
            value=data[1].lstrip(' ').rstrip(' ')
            parameters[parameter]=value
    #--
    return parameters


#Quantile
#input: list of values (floats/ints), p in [0,1)
#output: p-Quantile, 'n/a' if not applicable
def quantile(values,p):
    quantile='n/a'
    if values!=[]:
        dummy_values=[]
        for value in values:
            if value!='n/a':
                dummy_values.append(value)
        #--
        dummy_values.sort()
        num_values=len(dummy_values)
        index=int(float(num_values)*p)
        quantile=dummy_values[index]
    #-
    return quantile


#Median
#input: list of values (floats/ints)
#output: median, 'n/a' if not applicable
def median(values):
    median=quantile(values,0.5)
    return median


#generates a new file containing the given string
#input: string, filename
#output: -
def newFile(string,filename):
    file=open(filename,"w")
    file.write(string)
    file.close()

    
#executes the given program using the shell, writing output to stdout
#input: command (string)
#output: -
def run(args):
    process=Popen(args,shell=True)
    status=os.waitpid(process.pid, 0)
    return status[1]


#executes the given program using the shell, returning the output
#input: command (string)
#output: string
def run_return(args):
    process=Popen(args,shell=True,stdout=PIPE)
    return process.communicate()[0]


#wrapper function for running R scripts
#input: string of R commands
#output: stdout from R
def runR(R_commands):
    outfile_name='R_commands.dat'
    newFile(R_commands,outfile_name)
    run('R CMD BATCH '+outfile_name)
    infile_name=outfile_name+'.Rout'
    infile=open(infile_name)
    output=infile.read()
    return output


#Box-and-Whisker plot with median bars
#input: map x (float/int) -> list of y (float/int)
#output: list of (x,min,1st_quartile,median,3rd_quartile,max) tuples ready for gnuplot
def boxAndWhisker(data):
    tuples=[]
    for x in data:
        y_values=data[x]
        if y_values!=[]:
            minimum=min(y_values)
            maximum=max(y_values)
            median_value=median(y_values)
            first_quartile=quantile(y_values,0.25)
            third_quartile=quantile(y_values,0.75)
            tuples.append([x,minimum,first_quartile,median_value,third_quartile,maximum])
    #--
    return tuples


#cumulative probabilities
#input: list of values (float/int)), sort reverse (boolean)
#output: map: value -> cumulative probability
def cumulativeProbabilitiesMap(values,sort_reverse):
    cumulative_probabilities_map={}
    num_values=len(values)
    values_map={}
    for value in values:
        if value not in values_map:
            values_map[value]=1
        else:
            values_map[value]+=1
    #--
    sorted_values=sorted(values_map.keys(),reverse=sort_reverse)
    count=0
    for value in sorted_values:
        count+=values_map[value]
        cumulative_probabilities_map[value]=count/float(num_values)
    #-
    return cumulative_probabilities_map


#Histogram
#input: list of values (floats/ints), number of bins (int)
#output: list of pairs (bin,abundance) ready for gnuplot
def histogram(values,num_bins):
    cleaned_values=[]
    for value in values:
        if value!='n/a':
            cleaned_values.append(value)
    #--
    bins_list=[]
    maximum=max(cleaned_values)
    minimum=min(cleaned_values)
    #check if minimum unequal to maximum
    if maximum!=minimum:
        bin_size=float(maximum-minimum)/num_bins
        #init bins
        bins={}
        for j in range(num_bins):
            bins[j]=0
        #-
        #fill bins
        for value in cleaned_values:
            if value==maximum:
                bin_number=num_bins-1
            else:
                bin_number=math.floor(float(value-minimum)/float(bin_size))
            #-
            bins[bin_number]+=1
        #-
        #convert map of bin_numbers into list of bin representatives
        for bin_number in bins:
            bin=(bin_number*bin_size)+minimum
            bins_list.append((bin,bins[bin_number]))
    #--
    else:
        #create only one bin
        bins_list.append((minimum,len(values)))
    #-
    return bins_list


#Raw histogram
#input: list of values (floats/ints)
#output: list of pairs (bin,abundance) ready for gnuplot
def rawHistogram(values):
    cleaned_values=[]
    for value in values:
        if value!='n/a':
            cleaned_values.append(value)
    #--
    bins_list=[]
    bins={}
    #fill bins
    for value in cleaned_values:
        if value not in bins:
            bins[value]=1
        else:
            bins[value]+=1
    #--
    #convert map of bin_numbers into list of bin representatives
    for bin in bins:
        bins_list.append([bin,bins[bin]])
    #-
    return bins_list


#converts raw values into ranks for rank correlation coefficients
#input: list of values (int/float)
#output: map: value -> rank 
def getRanks(values):
    ranks={}
    sorted_values=sorted(values)
    for i in range(len(sorted_values)):
        value=sorted_values[i]
        if value not in ranks:
            ranks[value]=i+1
    #--
    return ranks


#Goodman and Kruskal's gamma correlation coefficient
#input: 2 lists of ranks (ints) of same length with corresponding entries
#output: Gamma correlation coefficient (rank correlation ignoring ties)
def gamma(ranks_list1,ranks_list2):
    num_concordant_pairs=0
    num_discordant_pairs=0
    num_tied_x=0
    num_tied_y=0
    num_tied_xy=0
    num_items=len(ranks_list1)
    for i in range(num_items):
        rank_1=ranks_list1[i]
        rank_2=ranks_list2[i]
        for j in range(i+1,num_items):
            diff1=ranks_list1[j]-rank_1
            diff2=ranks_list2[j]-rank_2
            if (diff1>0 and diff2>0) or (diff1<0 and diff2<0):
                num_concordant_pairs+=1
            elif (diff1>0 and diff2<0) or (diff1<0 and diff2>0):
                num_discordant_pairs+=1
            elif diff1==0 and diff2==0:
                num_tied_xy+=1
            elif diff1==0:
                num_tied_x+=1
            elif diff2==0:
                num_tied_y+=1
    #---
    try:
        gamma_corr_coeff=float(num_concordant_pairs-num_discordant_pairs)/float(num_concordant_pairs+num_discordant_pairs)
    except:
        gamma_corr_coeff='n/a'
    #-
    return [num_tied_x,num_tied_y,num_tied_xy,gamma_corr_coeff]


#Goodman and Kruskal's gamma correlation coefficient wrapper
#input: 2 lists of values of same length with corresponding entries
#output: Gamma correlation coefficient (rank correlation ignoring ties)
def gammaCC(values_list1,values_list2):
    ranks1=getRanks(values_list1)
    ranks_list1=[]
    for value in values_list1:
        rank=ranks1[value]
        ranks_list1.append(rank)
    #-
    ranks2=getRanks(values_list2)
    ranks_list2=[]
    for value in values_list2:
        rank=ranks2[value]
        ranks_list2.append(rank)
    #-
    gcc=round(gamma(ranks_list1,ranks_list2)[3],2)
    return gcc


#Change saturation of hexadecimal color
#input: hexadecimal color, saturation adjustment
#output: saturation-adjusted hexadecimal color string
def saturateHexColor(hexcolor, adjustment = 1.0):
    assert(adjustment >= 0 and len(hexcolor) >= 1)
    prefix = ""
    if hexcolor[0] == '#':
        hexcolor = hexcolor[1:]
        prefix = "#"
    #-
    assert(len(hexcolor) == 6)
    if adjustment == 1.0:
        return "%s%s" % (prefix, hexcolor)
    else:
        hsvColor = list(colorsys.rgb_to_hsv(int(hexcolor[0:2], 16)/255.0, int(hexcolor[2:4], 16)/255.0, int(hexcolor[4:6], 16)/255.0))
        hsvColor[1] = min(1.0, hsvColor[1] * adjustment)
        rgbColor = [min(255, 255 * v) for v in colorsys.hsv_to_rgb(hsvColor[0], hsvColor[1], hsvColor[2])]
        return "%s%.2x%.2x%.2x" % (prefix, rgbColor[0], rgbColor[1], rgbColor[2]) 


#Get list of hexadecimal colors
#input: number of colors needed (int), start value (optional), saturation adjustment (float; optional)
#output: list of hexadecimal color strings
def gnuplotColorWheel(n, start = 15, saturation_adjustment = None):
    hues = range(start, start + 360, 360/n)
    rgbcolors = ['%x%x%x' % (255 * hlscol[0], 255 * hlscol[1], 255 * hlscol[2]) for hlscol in [colorsys.hls_to_rgb(float(h % 360) / 360.0, 0.65, 1.00) for h in hues]]
    if saturation_adjustment:
        return [saturateHexColor(rgbcol, saturation_adjustment) for rgbcol in rgbcolors]
    else:
        return rgbcolors
