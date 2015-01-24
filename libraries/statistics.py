# This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

import math
import sys

#get bin size
#input: minimum (float/int), maximum (float/int), number of bins (int)
#output: bin size (float)
def getBinSize(minimum,maximum,num_bins):
    bin_size=float(maximum-minimum)/(num_bins-0.000001)#subtract epsilon to get the maximum also into a bin
    return bin_size


#init histogram
#input: minimum (float/int), maximum (float/int), bin size (float), num bins (int)
#output: map of representative value -> 0
def initHistogram(minimum,maximum,bin_size,num_bins):
    histogram={}
    if maximum!=minimum:
        for i in range(num_bins):
            binned_value=(i*bin_size)+minimum
            histogram[binned_value]=0
    #--
    else:
        histogram[minimum]=0
    #-
    return histogram


#converts value to representative value
#input: value (float/int), minimum (float/int), maximum (float/int), bin size (float)
#output: representative value (float/int)
def valueToBinnedValue(value,minimum,maximum,bin_size):
    #check if minimum unequal to maximum
    if maximum!=minimum:
        binned_value=math.floor((value-minimum)/float(bin_size))*bin_size+minimum
    else:
        binned_value=minimum
    #-
    return binned_value


#init log histogram
#input: minimum (float/int) !=0, maximum (float/int) !=0, log bin size (float)
#output: map of representative value (log bin size = number of potencies of 10 to merge) -> 0
def initLogHistogram(minimum,maximum,log_bin_size):
    histogram={}
    min_bin=int(math.ceil(math.floor(math.log10(minimum))/(-log_bin_size))*(-log_bin_size))
    max_bin=int(math.ceil(math.floor(math.log10(maximum))/(-log_bin_size))*(-log_bin_size))
    bin=min_bin
    while bin<=max_bin:
        representative_value=math.pow(10,bin)
        histogram[representative_value]=0
        bin+=log_bin_size
    #-
    histogram[0]=0
    return histogram


#converts log value to representative value
#input: log value (float/int), minimum (float/int), maximum (float/int), bin size (float)
#output: representative value (float/int)
def logValueToBinnedValue(log_value,log_bin_size):
    if log_value>0:
        bin=int(math.ceil(math.floor(math.log10(log_value))/(-log_bin_size))*(-log_bin_size))
        binned_value=math.pow(10,bin)
    else:
        binned_value=0
    #-
    return binned_value


#standard binning of values (stable bin size)
#input: value (float/int), bin_size (float)
#output: bin (float/int)
def getBin(value,bin_size):
    bin=math.floor(value/float(bin_size))*bin_size
    return bin


#standard binning of values (stable bin size)
#input: value (float/int), list of values (floats/ints), number of bins (int)
#output: bin (float/int)
def getBin2(value,values,num_bins):
    cleaned_values=[]
    for value in values:
        if value!='n/a':
            cleaned_values.append(value)
    #--
    maximum=max(cleaned_values)
    minimum=min(cleaned_values)
    #check if minimum unequal to maximum
    if maximum!=minimum:
        bin_size=float(maximum-minimum)/(num_bins-0.000001)#subtract epsilon to get the maximum also into a bin
        bin=math.floor(value/float(bin_size))*bin_size
    else:
        bin=minimum
    #-
    return bin


#raw to standard histogram
#input: raw value histogram (map: float -> int), min, max, number of bins
#output: map: bin (float) -> abundance (int), mapping of original to binned values
def rawToStandardHistogram(raw_value_histogram,minimum,maximum,num_bins):
    histogram={}
    mapping={}
    bins={}
    #init standard bins
    bins[0]=0
    #check if minimum unequal to maximum
    if maximum!=minimum:
        for i in range(1,num_bins):
            bins[i]=0
    #--
    bin_size=float(maximum-minimum)/(num_bins-0.000001)#subtract epsilon to get the maximum also into a bin
    #convert raw values into bins
    for value in raw_value_histogram:
        bin=0
        #check if minimum unequal to maximum
        if maximum!=minimum:
            bin=math.floor(float(value-minimum)/float(bin_size))
        #-
        bins[bin]+=raw_value_histogram[value]
        representative_value=(bin*bin_size)+minimum
        mapping[value]=representative_value
    #-
    #generate histogram of representative values
    for bin in bins:
        representative_value=(bin*bin_size)+minimum
        histogram[representative_value]=bins[bin]
    #-
    return (histogram,mapping)


#raw to log histogram (special treatment of zero)
#input: raw value histogram (map: float -> int), min!=0, max!=0, bin size 
#output: map: bin (float) -> (abundance (int) with a bin for every bin_size potency of 10, mapping of original to binned values)
def rawToLogHistogram(raw_value_histogram,minimum,maximum,bin_size):
    histogram={}
    mapping={}
    bins={}
    #init logarithmic bins
    min_bin=int(math.ceil(math.floor(math.log10(minimum))/(-bin_size))*(-bin_size))
    max_bin=int(math.ceil(math.floor(math.log10(maximum))/(-bin_size))*(-bin_size))
    bin=min_bin
    while bin<=max_bin:
        bins[bin]=0
        bin+=bin_size
    #-
    #convert raw values into bins
    for value in raw_value_histogram:
        bin='n/a'
        if value>0:
            bin=int(math.ceil(math.floor(math.log10(value))/(-bin_size))*(-bin_size))
        #-
        if bin=='n/a':
            bins[bin]=raw_value_histogram[value]
            representative_value=0
        else:
            bins[bin]+=raw_value_histogram[value]
            representative_value=math.pow(10,bin)
        #-
        mapping[value]=representative_value            
    #--
    #generate histogram of representative values
    for bin in bins:
        if bin=='n/a':
            representative_value=0
        else:
            representative_value=math.pow(10,bin)
        #-
        histogram[representative_value]=bins[bin]
    #-
    return (histogram,mapping)


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
        bin_size=float(maximum-minimum)/(num_bins-0.000001)#subtract small epsilon value to get the maximum also into a bin
        #init bins
        bins={}
        for j in range(num_bins):
            bins[j]=0
        #-
        #fill bins
        for value in cleaned_values:
            bin_number=math.floor(float(value-minimum)/float(bin_size))
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


#Histogram2
#input: list of values (floats/ints), number of bins (int)
#output: list of pairs (bin,abundance) ready for gnuplot
def histogram2(values,num_bins,epsilon):
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
        bin_size=float(maximum-minimum)/(num_bins-epsilon)#subtract epsilon to get the maximum also into a bin
        #init bins
        bins={}
        for j in range(num_bins):
            bins[j]=0
        #-
        #fill bins
        for value in cleaned_values:
            bin_number=math.floor(float(value-minimum)/float(bin_size))
            if bin_number==num_bins:
                bins[bin_number-1]+=1
            else:
                bins[bin_number]+=1
        #--
        #convert map of bin_numbers into list of bin representatives
        for bin_number in bins:
            bin=(bin_number*bin_size)+minimum
            bins_list.append([bin,bins[bin_number]])
    #--
    else:
        #create only one bin
        bins_list.append([minimum,len(values)])
    #-
    return bins_list


#log histogram
#input: list of values (floats/ints)
#output: list of pairs (bin,abundance) with a bin for every potence of 10
def logHistogram(values):
    bins_list=[]
    bins={}
    #fill bins
    for value in values:
        if value==0:
            bin_number=0
        else:
            bin_number=int(math.floor(math.log10(value)))
        #-
        if bin_number not in bins:
            bins[bin_number]=1
        else:
            bins[bin_number]+=1
    #--
    #convert map of bin numbers into list of bin representatives
    for bin_number in bins:
        if bin_number==0:
            bin=0
        else:
            bin=math.pow(10,bin_number)
        #-
        bins_list.append([bin,bins[bin_number]])
    #-
    return bins_list


#Average
#input: list of values
#output: average value, 'n/a' if not applicable
def average(values):
    average='n/a'
    if len(values)!=0:
        average=0
        for value in values:
            average+=float(value)
        #-
        average/=len(values)
    #-
    return average


#Average and standard deviation
#input: list of values
#output: tupel (average value,standard deviation) ('n/a' if not applicable)
def averageAndStddev(values):
    average='n/a'
    stddev='n/a'
    num_values=len(values)
    if num_values!=0:
        average=0
        for value in values:
            average+=float(value)
        #-
        average/=float(num_values)
        variance=0
        for value in values:
            dev=float(value)-average
            variance+=float(dev*dev)
        #-
        variance/=float(num_values)
        stddev=math.sqrt(variance)
    #-
    return (average,stddev)


#Average and standard error
#input: list of values
#output: tupel (average value,standard error) ('n/a' if not applicable)
def averageAndStderr(values):
    average='n/a'
    stderr='n/a'
    num_values=len(values)
    if num_values!=0:
        average=0
        for value in values:
            average+=float(value)
        #-
        average/=float(num_values)
        variance=0
        for value in values:
            dev=float(value)-average
            variance+=float(dev*dev)
        #-
        if num_values>1:
            variance/=float(num_values-1)#Bessel's correction
        else:
            variance/=float(num_values)
        #-
        stddev=math.sqrt(variance)
        stderr=stddev/math.sqrt(num_values)
    #-
    return (average,stderr)    


#Average and standard deviation
#input: map of value (float/int) -> occurrence (int)
#output: tupel (num_values,average value,standard deviation) ('n/a' if not applicable)
def averageAndStddev2(map):
    average='n/a'
    stddev='n/a'
    num_values=sum(map.values())
    if num_values!=0:
        average=0
        for value in map:
            average+=map[value]*value
        #-
        average/=float(num_values)
        variance=0
        for value in map:
            dev=value-average
            variance+=map[value]*(dev*dev)
        #-
        variance/=float(num_values)
        stddev=math.sqrt(variance)
    #-
    return (num_values,average,stddev)
    

#Enrichment factor
#input: value1, value2
#output: enrichment factor, 'n/a' if not applicable
def enrichmentFactor(value1,value2):
    enrichment_factor='n/a'
    if float(value2)!=0:
        enrichment_factor=float(value1)/float(value2)
    #-
    return enrichment_factor


#Enrichment
#input: value1, value2
#output: enrichment value, 'n/a' if not applicable
def enrichment(value1,value2):
    enrichment='n/a'
    enrichment_factor=enrichmentFactor(value1,value2)
    if enrichment_factor!=0 and enrichment_factor!='n/a':
        #take the log_2 to get a symmetric range
        enrichment=math.log(enrichment_factor,2)
    #-
    return enrichment


#Calculation of enrichment
#input: list of values, list/map of target values, list/map of reference values
#output: float
def calculateEnrichment(values,target_set,reference_set):
    target_set_overlap_size=0
    reference_set_overlap_size=0
    for value in values:
        if value in target_set:
            target_set_overlap_size+=1
        #-
        if value in reference_set:
            reference_set_overlap_size+=1
    #--
    target_set_size=len(target_set)
    reference_set_size=len(reference_set)
    fraction_target_set_covered=target_set_overlap_size/float(target_set_size)
    fraction_reference_set_covered=reference_set_overlap_size/float(reference_set_size)
    return enrichment(fraction_target_set_covered,fraction_reference_set_covered)


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


#Median (upper median)
#input: list of values (floats/ints)
#output: median, 'n/a' if not applicable
def median(values):
    median=quantile(values,0.5)
    return median


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


#Tukey Box-and-Whisker plot with median bars, whiskers denoting +/- 1.5xIQR, outliers are separate
#input: map x (float/int) -> list of y (float/int)
#output: map of x (float/int) -> ([x,lower,1st_quartile,median,3rd_quartile,upper],[outliers]) tuples ready for gnuplot
def tukeyBoxAndWhisker(data):
    from collections import namedtuple

    Stats = namedtuple( 'Stats', [
                'x',
                'lower_whisker',
                'first_quartile',
                'median',
                'third_quartile',
                'upper_whisker'])

    results={}
    for x in data:
        y_values=data[x]
        if y_values!=[]:
            median_value=median(y_values)
            first_quartile=quantile(y_values,0.25)
            third_quartile=quantile(y_values,0.75)
            interquartile_range=third_quartile-first_quartile
            lower=first_quartile
            upper=third_quartile
            for value in y_values:
                if value<lower and value>=first_quartile-1.5*interquartile_range:
                    lower=value
                elif value>upper and value<=third_quartile+1.5*interquartile_range:
                    upper=value
            #--
            outliers=[]
            for value in y_values:
                if value<lower or value>upper:
                    outliers.append(value)
            #--
            results[x] = Stats(x, lower, first_quartile, median_value, third_quartile, upper), outliers
    #--
    return results


#p-value
#input: value, list of background values
#output: probability of getting this value by chance (n/a if not applicable)
def pValue(value,background):
    p_value='n/a'
    if len(background)==0 or value=='n/a':
        return p_value
    #-
    count=0
    value=float(value)
    if value>=0:
        for item in background:
            if item!='n/a' and float(item)>=value:
                count+=1
    #---
    else:
        for item in background:
            if item!='n/a' and float(item)<=value:
                count+=1
    #---
    if count==0:
        p_value='<'+str(1/float(len(background)))
    else:
        p_value=float(count)/float(len(background))
    #-
    return p_value


#interquartile range
#input: list of values (floats/ints)
#output: Q3-Q1
def interquartileRange(values):
    iqrange='n/a'
    third_quartile=quantile(values,0.75)
    first_quartile=quantile(values,0.25)
    if third_quartile!='n/a' and first_quartile!='n/a':
        iqrange=third_quartile-first_quartile
    #-
    return iqrange


#Z-score standardization using median and interquartile range to avoid bias from outliers
#input: list of values (floats/ints)
#output: standardized list of values
def zScoreStandardization(values):
    standardized_values=[]
    median_value=median(values)
    interquartile_range=interquartileRange(values)
    iq_range_norm=1.34898
    if median_value!='n/a' and interquartile_range!='n/a':
        for value in values:
            if value!='n/a':
                standardized_values.append((value-median_value)/float(interquartile_range/iq_range_norm))
            else:
                standardized_values.append('n/a')
    #--
    return standardized_values


#euklidian distance
#input: 2 lists of values (floats/ints) of same length
#output: euklidian distance between the two lists
def euklidianDistance(profile1,profile2):
    distance=0
    num_conditions=len(profile1)
    if num_conditions!=len(profile2) or 'n/a' in profile1 or 'n/a' in profile2:
        return 'n/a'
    #-
    for i in range(num_conditions):
        x=float(profile1[i])
        y=float(profile2[i])
        dif=x-y
        distance+=dif*dif
    #-
    distance=math.sqrt(distance)
    return distance


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


#Goodman and Kruskal's gamma correlation coefficient wrapper 2
#input: set of relevant items, 2 maps of item -> value (relevant items have to be part of the maps),outfile_name
#output: Gamma correlation coefficient (rank correlation ignoring ties) rounded to two decimals
def gcc(relevant_items,value_map1,value_map2,outfile_name=None):
    values_list1=[]
    values_list2=[]
    for item in relevant_items:
        values_list1.append(value_map1[item])
        values_list2.append(value_map2[item])
    #-
    num_values=len(values_list1)
    ranks1=getRanks(values_list1)
    ranks2=getRanks(values_list2)
    ranks_list1=[]
    ranks_list2=[]
    for item in relevant_items:
        value1=value_map1[item]
        value2=value_map2[item]
        ranks_list1.append(ranks1[value1])
        ranks_list2.append(ranks2[value2])
    #-
    if outfile_name!=None:
        outfile=open(outfile_name,"w")
        for i in range(num_values):
            outfile.write(str(values_list1[i])+'\t'+str(values_list2[i])+'\n')
        #-
        outfile.close()
    #-
    gcc=round(gamma(ranks_list1,ranks_list2)[3],2)
    return gcc


#Kendall's tau_a correlation coefficient
#input: 2 lists of ranks (ints) of same length
#output: tau_a correlation coefficient (rank correlation)
def tau_a(ranks1,ranks2):
    num_concordant_pairs=0
    num_discordant_pairs=0
    num_items=len(ranks1)
    for i in range(num_items):
        for j in range(i+1,num_items):
            diff1=ranks1[j]-ranks1[i]
            diff2=ranks2[j]-ranks2[i]
            if (diff1>0 and diff2>0) or (diff1<0 and diff2<0):
                num_concordant_pairs+=1
            elif (diff1>0 and diff2<0) or (diff1<0 and diff2>0):
                num_discordant_pairs+=1
    #---
    total_num_pairs=float(num_items*(num_items-1))/float(2)
    return float(num_concordant_pairs-num_discordant_pairs)/float(total_num_pairs)


#Kendall's tau_b correlation coefficient
#input: 2 lists of ranks (ints) of same length
#output: tau_b correlation coefficient (rank correlation)
def tau_b(ranks1,ranks2):
    num_concordant_pairs=0
    num_discordant_pairs=0
    num_tied_x=0
    num_tied_y=0
    num_tied_xy=0
    num_items=len(ranks1)
    for i in range(num_items):
        for j in range(i+1,num_items):
            diff1=ranks1[j]-ranks1[i]
            diff2=ranks2[j]-ranks2[i]
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
    return float(num_concordant_pairs-num_discordant_pairs)/float(math.sqrt(float((num_concordant_pairs+num_discordant_pairs+num_tied_x)*(num_concordant_pairs+num_discordant_pairs+num_tied_y))))


#Pearson's R (sample) correlation coefficient
#input: 2 lists of values (x,y) of same length with corresponding entries
#output: Pearson correlation coefficient (requires normally distributed data)
def pearsonCC(x,y):
    assert len(x) == len(y)
    assert len(x) > 0
    avg_x=average(x)
    avg_y=average(y)
    diffprod_sum=0
    xdiff2_sum=0
    ydiff2_sum=0
    for i in range(len(x)):
        xdiff=x[i]-avg_x
        ydiff=y[i]-avg_y
        diffprod_sum+=xdiff*ydiff
        xdiff2_sum+=xdiff*xdiff
        ydiff2_sum+=ydiff*ydiff
    #-
    return round(diffprod_sum/math.sqrt(xdiff2_sum*ydiff2_sum),2)


#false-positive rate
#input: FP,TN (ints)
#output: float
def false_positive_rate(FP,TN):
    false_positive_rate='n/a'
    if FP+TN!=0:
        false_positive_rate=FP/float(FP+TN)
    #-
    return false_positive_rate


#true-positive rate
#input: TP,FN (ints)
#output: float
def true_positive_rate(TP,FN):
    true_positive_rate='n/a'
    if TP+FN!=0:
        true_positive_rate=TP/float(TP+FN)
    #-
    return true_positive_rate


#false-negative rate
#input: FN,TP (ints)
#output: float
def false_negative_rate(FN,TP):
    false_negative_rate='n/a'
    if FN+TP!=0:
        false_negative_rate=FN/float(FN+TP)
    #-
    return false_negative_rate


#sensitivity: TP/(TP+FN)
#input: TP,FN (ints)
#output: float
def sensitivity(TP,FN):
    sensitivity='n/a'
    if TP+FN!=0:
        sensitivity=TP/float(TP+FN)
    #-
    return sensitivity


#specificity: TN/(TN+FP)
#input: TN,FP (ints)
#output: float
def specificity(TN,FP):
    specificity='n/a'
    if TN+FP!=0:
        specificity=TN/float(TN+FP)
    #-
    return specificity


#false-discovery rate (informal): FP/(FP+TP)
#input: FP,TP (ints)
#output: float
def falseDiscoveryRate(FP,TP):
    false_discovery_rate='n/a'
    if FP+TP!=0:
        false_discovery_rate=FP/float(FP+TP)
    #-
    return false_discovery_rate


#positive predictive value: TP/(TP+FP)
#input: TP,FP (ints)
#output: float
def positivePredictiveValue(TP,FP):
    positive_predictive_value='n/a'
    if TP+FP!=0:
        positive_predictive_value=TP/float(TP+FP)
    #-
    return positive_predictive_value


#precision: TP/(TP+FP)
#input: TP,FP (ints)
#output: float
def precision(TP,FP):
    precision='n/a'
    if TP+FP!=0:
        precision=TP/float(TP+FP)
    #-
    return precision


#recall: TP/(TP+FN)
#input: TP,FN (ints)
#output: float
def recall(TP,FN):
    recall='n/a'
    if TP+FN!=0:
        recall=TP/float(TP+FN)
    #-
    return recall


#F-measure: harmonic mean of precision and recall
#input: precision,recall
#output: float
def fMeasure(precision,recall):
    f_measure='n/a'
    if precision!='n/a' and recall!='n/a' and precision+recall!=0:
        f_measure=2*precision*recall/(precision+recall)
    #-
    return f_measure


#negative predictive value: TN/(TN+FN)
#input: TN,FN (ints)
#output: float
def negativePredictiveValue(TN,FN):
    negative_predictive_value='n/a'
    if TN+FN!=0:
        negative_predictive_value=TN/float(TN+FN)
    #-
    return negative_predictive_value


#accuracy: (TP+TN)/(TP+TN+FP+FN)
#input: TP,TN,FP,FN (ints)
#output: float
def accuracy(TP,TN,FP,FN):
    accuracy='n/a'
    if TP+TN+FP+FN!=0:
        accuracy=(TP+TN)/float(TP+TN+FP+FN)
    #-
    return accuracy


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


#cumulative probabilities
#input: map of value (float/int) -> occurrence (int), sort reverse (boolean)
#output: map: value -> cumulative probability
def cumulativeProbabilitiesMap2(values_map,sort_reverse):
    cumulative_probabilities_map={}
    num_values=sum(values_map.values())
    sorted_values=sorted(values_map.keys(),reverse=sort_reverse)
    count=0
    for value in sorted_values:
        count+=values_map[value]
        cumulative_probabilities_map[value]=count/float(num_values)
    #-
    return cumulative_probabilities_map


#log-odds ratio
#input: fraction of test-items fulfilling the criterion, fraction of control-items fulfilling the criterion
#output: float, 'n/a' if not applicable
def logOddsRatio(fraction_test_items,fraction_control_items):
    log_odds_ratio='n/a'
    if fraction_control_items!=0 and (1-fraction_test_items)!=0:
        odds_ratio=(fraction_test_items*(1-fraction_control_items)/float((1-fraction_test_items)*fraction_control_items))
        if odds_ratio!=0:
            log_odds_ratio=math.log(odds_ratio)
    #--
    return log_odds_ratio


#node degree
#input: node reference id, node-based network map
#output: int, 'n/a' if the node id cannot be found
def nodeDegree(node_id,network):
    node_degree='n/a'
    if node_id in network:
        node_degree=len(network[node_id])
    #-
    return node_degree


#local clustering coefficient for a given node (for undirected graphs, ignoring self-edges)
#input: node reference id, node-based network map
#output: float (fraction of all possible links that exist between the direct neighbors of the given node), 'n/a' if the node id cannot be found
def localClusteringCoefficient(node_id,network):
    clustering_coefficient='n/a'
    if node_id in network:
        interactors=network[node_id]
        #determine list of direct neighbors
        direct_neighbors=[]
        for interactor in interactors:
            if interactor!=node_id:
                direct_neighbors.append(interactor)
        #--
        num_direct_neighbors=len(direct_neighbors)
        num_possible_links=num_direct_neighbors*(num_direct_neighbors-1)/2
        num_existing_links=0
        for i in range(num_direct_neighbors):
            node1=direct_neighbors[i]
            for j in range(i+1,num_direct_neighbors):
                node2=direct_neighbors[j]
                if node2 in network[node1]:
                    num_existing_links+=1
        #---
        if num_possible_links!=0:
            clustering_coefficient=num_existing_links/float(num_possible_links)
    #--
    return clustering_coefficient


#Chi-square test (ignores values which occur less than 5 times)
#input: list of ints (values1), list of ints (values2)
#output: p-value (float)
def chiSquareTest(values1,values2):
    p_value='n/a (less than 2 values per sample)'
    import rpy
    min_value=min(min(values1),min(values2))
    max_value=max(max(values1),max(values2))
    all_values=[]
    for i in range(min_value,max_value+1):
        x=values1.count(i)
        y=values2.count(i)
        if x>5 and y>5:
            all_values.append(x)
            all_values.append(y)
    #--
    #test if enough values for 2x2 contingency table
    if len(all_values)>=4:
        matrix=rpy.r.matrix(all_values,nr=2)
        p_value=rpy.r.chisq_test(matrix)['p.value']
    #-
    return p_value
    

#Chi-square test (ignores values which occur less than 5 times)
#input: two maps of value (float/int) -> occurrence (int)
#output: p-value (float)
def chiSquareTest2(map1,map2):
    p_value='n/a (less than 2 values per sample)'
    import rpy
    all_keys=map1.keys()
    all_keys.extend(map2.keys())
    all_keys_set=set(all_keys)
    sorted_all_keys_list=sorted(all_keys_set)
    all_values=[]
    for key in sorted_all_keys_list:
        x=0
        if key in map1:
            x=map1[key]
        #-
        y=0
        if key in map2:
            y=map2[key]
        #-
        if x>5 and y>5:
            all_values.append(x)
            all_values.append(y)
    #--
    #test if enough values for 2x2 contingency table
    if len(all_values)>=4:
        matrix=rpy.r.matrix(all_values,nr=2)
        p_value=rpy.r.chisq_test(matrix)['p.value']
    #-
    return p_value


#Wilcoxon rank sum test (Mann-Whitney U test)
#input: two maps of value (float/int) -> occurrence (int), one of {'two.sided','greater','less'}
#output: one or two-sided p-value (float), depending on the given alternative
def wilcoxonTest(map1,map2,alternative):
    p_value='n/a (less than 5 values per sample)'
    import rpy
    x=[]
    for key in map1:
        occurrence=map1[key]
        for i in range(occurrence):
            x.append(key)
    #--
    y=[]
    for key in map2:
        occurrence=map2[key]
        for i in range(occurrence):
            y.append(key)
    #--
    if alternative=='two.sided':
        print 'two-sided test'
    else:
        print 'one-sided test, alternative='+alternative
    #-
    if len(x)>=5 and len(y)>=5:
        p_value=rpy.r.wilcox_test(x,y,alternative)['p.value']
    #-
    return p_value


#Wilcoxon rank sum test (Mann-Whitney U test)
#input: two lists of values (float/int), one of {'two.sided','greater','less'}
#output: one or two-sided p-value (float), depending on the given alternative
def wilcoxonTest2(list1,list2,alternative):
    p_value='n/a (less than 5 values per sample)'
    import rpy
    x=[]
    x.extend(list1)
    y=[]
    y.extend(list2)
    if alternative=='two.sided':
        print 'two-sided test'
    else:
        print 'one-sided test, alternative='+alternative
    #-
    if len(x)>=5 and len(y)>=5:
        p_value=rpy.r.wilcox_test(x,y,alternative)['p.value']
    #-
    return p_value


#Kolmogorov-Smirnov test (KS test)
#input: two lists of values (float/int), one of {'two.sided','greater','less'}
#output: one or two-sided p-value (float), depending on the given alternative
def ksTest(list1,list2,alternative):
    p_value='n/a'
    import rpy
    x=[]
    x.extend(list1)
    y=[]
    y.extend(list2)
    if alternative=='two.sided':
        print 'two-sided test'
    else:
        print 'one-sided test, alternative='+alternative
    #-
    #if len(x)>=5 and len(y)>=5:
    p_value=rpy.r.ks_test(x,y,alternative)['p.value']
    #-
    return p_value


#t-test
#input: two lists of values (float/int), one of {'two.sided','greater','less'}
#output: one or two-sided p-value (float), depending on the given alternative
def tTest(list1,list2,alternative):
    p_value='n/a'
    import rpy
    x=[]
    x.extend(list1)
    y=[]
    y.extend(list2)
    if alternative=='two.sided':
        print 'two-sided test'
    else:
        print 'one-sided test, alternative='+alternative
    #-
    #if len(x)>=5 and len(y)>=5:
    p_value=rpy.r.t_test(x,y,alternative)['p.value']
    #-
    return p_value


#Fisher's exact test
#input: two lists of values (float/int), one of {'two.sided','greater','less'}, one of {True,False} (optional)
#output: one or two-sided p-value (float), depending on the given alternative
def fishersExactTest(list1,list2,alternative,simulate=False,printout=True):
    p_value='n/a'
    import rpy
    values=[]
    for i in range(len(list1)):
        values.append(list1[i])
        values.append(list2[i])
    #-
    matrix=rpy.r.matrix(values,nrow=2)
    p_value_info=rpy.r.fisher_test(matrix,workspace=5.2e8,alternative=alternative,simulate_p_value=simulate)
    if printout:
        if alternative=='two.sided':
            print 'two-sided test'
        else:
            print 'one-sided test, alternative='+alternative
        #-
        print matrix
        print p_value_info['alternative']
    #--
    p_value=p_value_info['p.value']
    return p_value


#Fisher's exact test for nxm contingency tables
#input: list of lists of values (float/int), one of {'two.sided','greater','less'}, one of {True,False} (optional)
#output: one or two-sided p-value (float), depending on the given alternative
def fishersExactTest2(list_of_lists,alternative,simulate=False):
    p_value='n/a'
    import rpy
    num_lists=len(list_of_lists)
    if alternative=='two.sided':
        print 'two-sided test'
    else:
        print 'one-sided test, alternative='+alternative
    #-
    values=[]
    list1=list_of_lists[0]
    for i in range(len(list1)):
        for j in range(num_lists):
            values.append(list_of_lists[j][i])
    #--
    matrix=rpy.r.matrix(rpy.r.c(values),nr=num_lists)
    print matrix
    p_value_info=rpy.r.fisher_test.lcall((('x',matrix),('workspace',5.2e8),('alternative',alternative),('simulate',simulate)))
    print p_value_info['alternative']
    p_value=p_value_info['p.value']
    return p_value


#Isotone regression
#input: map of key (float/int) -> value (float/int), regression type (i.e. one of 'isotonic' and 'antitonic'
#output: map of key -> fitted value
def isotoneRegression(map,regression_type):
    import rpy
    print 'loading R library "centered isotone regression (CIR)"...'
    rpy.r.library('cir')
    fitted_data={}
    x=[]
    y=[]
    sorted_keys=sorted(map.keys())
    for key in sorted_keys:
        x.append(key)
        y.append(map[key])
    #-
    if regression_type=='antitonic':
        print 'regression type: antitonic'
        monotone_regression=rpy.r.pava(y,dec=True)
    else:
        print 'regression type: isotonic'
        monotone_regression=rpy.r.pava(y,dec=False)
    #-
    x_values=x
    fitted_y_values=monotone_regression
    for i in range(len(x_values)):
        key=x_values[i]
        fitted_value=fitted_y_values[i]
        fitted_data[key]=fitted_value
    #-
    return fitted_data


#Centered isotone regression (suited for data which is assumed to be strictly monotone; avoids flat stretches)
#input: map of key (float/int) -> value (float/int), regression type (i.e. one of 'isotonic' and 'antitonic'
#output: map of key -> fitted value
def centeredIsotoneRegression(map,regression_type):
    import rpy
    print 'loading R library "centered isotone regression (CIR)"...'
    rpy.r.library('cir')
    fitted_data={}
    x=[]
    y=[]
    sorted_keys=sorted(map.keys())
    for key in sorted_keys:
        x.append(key)
        y.append(map[key])
    #-
    if regression_type=='antitonic':
        print 'regression type: antitonic'
        monotone_regression=rpy.r.cir_pava(y,x,dec=True,full=True)
    else:
        print 'regression type: isotonic'
        monotone_regression=rpy.r.cir_pava(y,x,dec=False,full=True)
    #-
    x_values=monotone_regression['original.x']
    fitted_y_values=monotone_regression['output.y']
    for i in range(len(x_values)):
        key=x_values[i]
        fitted_value=fitted_y_values[i]
        fitted_data[key]=fitted_value
    #-
    return fitted_data


#Non-linear least-squares fit to Gaussian curve
#input: data map of x value -> corresponding y value, bool (print trace)
#output: (residual sum of square, map of parameter -> fitted value)
def nlsGaussianFit(data,print_trace):
    rss='n/a'
    parameters_map={}
    import rpy
    x_values=[]
    y_values=[]
    sorted_keys=sorted(data.keys())
    for item in sorted_keys:
        x_values.append(item)
        y_values.append(data[item])
    #-
    amplitude=max(y_values)
    mean=0
    stddev=1
    for i in range(len(x_values)):
        if y_values[i]==amplitude:
            mean=x_values[i]
            break
    #--
    data_frame=rpy.r.data_frame(x=x_values,y=y_values)
    start_frame=rpy.r.data_frame(a=amplitude,m=mean,s=stddev)
    #control_frame=rpy.r.nls_control(minFactor = 1/2048, maxiter = 1000,tol = 1e-04)
    nls_fit=rpy.r.nls('y ~ a/s*exp(-(x-m)**2/(2*s**2))',data=data_frame,start=start_frame,trace=print_trace, algorithm='port')
    nls_fit_results=nls_fit['m']
    if nls_fit_results!=None:
        rss=nls_fit_results['deviance']()
        parameters=nls_fit_results['getPars']()
        for item in parameters:
            parameters_map[item]=parameters[item]
    #--
    return (rss,parameters_map)
    

#Normal distribution
#input: x value, scaling factor, mean and standard deviation, computed by Gnuplot fit
#output: y value
def normalDistribution(x,a,mu,sigma):
    y_value='n/a'
    if sigma!=0:
        y_value=a/float(sigma)*math.exp(-(x-mu)**2/float(2*sigma**2))
    #-
    return y_value


#Monotone regression using the pool adjacent violators algorithm (pava)
#input: map of key (float/int) -> value (float/int), regression type (i.e. one of 'isotonic' and 'antitonic')
#output: map of key -> fitted value
def monotoneRegression(map,regression_type):
    fitted_data={}
    sorted_keys=sorted(map.keys())
    values=[]
    for key in sorted_keys:
        values.append(map[key])
    #-
    num_values=len(values)
    if regression_type=='antitonic':
        for i in range(num_values):
            values[i]*=-1
    #--
    levels=[]
    for i in range(num_values):
        levels.append(i)
    #-
    #print 'values '+str(values)
    #print 'levels '+str(levels)
    while True:
        #find adjacent violators
        no_violators=True
        adjacent_violators=[]
        for i in range(num_values-1):
            if values[i+1]-values[i]<0:
                #violator found
                adjacent_violators.append(True)
                no_violators=False
            else:
                adjacent_violators.append(False)
        #--
        #print 'adjacent violators '+str(adjacent_violators)
        if no_violators:
            break
        #-
        #pool first pair of violators
        for i in range(num_values-1):
            if adjacent_violators[i]:
                index_first_violator=i
                break
        #--
        #print 'index first violator '+str(index_first_violator)
        level1=levels[index_first_violator]
        level2=levels[index_first_violator+1]
        #print 'level1 '+str(level1)
        #print 'level2 '+str(level2)
        level_flags=[]
        for i in range(num_values):
            level_flags.append(levels[i]==level1 or levels[i]==level2)
        #-
        #print 'level-flags '+str(level_flags)
        adjusted_value=0
        count=0
        for i in range(num_values):
            if level_flags[i]:
                count+=1
                adjusted_value+=values[i]
        #--
        adjusted_value/=float(count)
        for i in range(num_values):
            if level_flags[i]:
                values[i]=adjusted_value
        #--
        #print 'values '+str(values)
        for i in range(num_values):
            if level_flags[i]:
                levels[i]=level1
        #--
        #print 'levels '+str(levels)
    #-
    #print
    #print 'final values '+str(values)
    if regression_type=='antitonic':
        for i in range(num_values):
            values[i]*=-1
    #--    
    for i in range(num_values):
        key=sorted_keys[i]
        fitted_value=values[i]
        fitted_data[key]=fitted_value
    #-
    return fitted_data
