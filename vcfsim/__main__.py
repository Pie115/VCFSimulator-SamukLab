import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
import argparse
from IPython.display import SVG, display
from .SimulatorClass import MyVcfSim
import warnings

def multiple_chrom(chromfilename = 'input.txt', seed = 1234, foldername = 'PixyFolder', percentmissing = 0, percentsitemissing = 0, outputfile = 'myvcftest', samp_num = 20, sample_names = None):
    
    f = open(chromfilename, "r")

    a = f.readline()
    
    chromlist = []
    ploidylist = []
    lengthlist = []
    nelist = []
    mulist = []

    #check if all lists are equal

    while(a):
        a = a.split()
        if(len(a) != 5):
            print('Error, sample file must only have 5 columns, chromosone, ploidy, sequence length, population size, mutation rate')
            sys.exit(0)
            
        if(not a[1].isdigit()):
            print('Error, ploidy must be a number')
            sys.exit(0)

        try:
            temp = int(eval(a[2]))
        except Exception as e:
            print("Error, sequence length must be an expression or number")
            sys.exit(0)

        try:
            temp = int(eval(a[3]))
        except Exception as e:
            print("Error, effective population must be an expression or number")
            sys.exit(0)
            
        try:
            temp = int(eval(a[4]))
        except Exception as e:
            print("Error, mutation rate must be an expression or number")
            sys.exit(0)
            

        chromlist.append(a[0]), ploidylist.append(int(a[1])), lengthlist.append(int(eval(a[2]))), nelist.append(int(eval(a[3]))),  
        mulist.append(eval(a[4]))
        
        a = f.readline()

    temp = len(chromlist) + len(ploidylist) + len(lengthlist) + len(nelist) + len(mulist)
    if(temp/5 != len(chromlist)):
        print('Error, all columns must be the same length')
        sys.exit(0)

    for i in range(len(chromlist)):
        sim = MyVcfSim(chromlist[i], lengthlist[i], ploidylist[i], nelist[i], mulist[i], percentmissing, percentsitemissing, seed, 
                       chromlist[i], samp_num, 'population.txt', 'vcf', sample_names)
        
        sim.simulate_vcfs()
    
    with open(outputfile, "w") as f:
        pass

    for i in range(len(chromlist)):
        if i == 0:
            with open(chromlist[i], 'r') as source:
                lines = source.readlines()
            with open(outputfile, 'a') as destination:
                destination.writelines(lines)
        else:
            with open(chromlist[i], 'r') as source:
                lines = source.readlines()
            with open(outputfile, 'a') as destination:
                destination.writelines(lines[6:])
        os.remove(chromlist[i])

def vcf_simulator(chrom = 1, amountofruns = 1, seed = 1234, foldername = 'PixyFolder', sitesize = 10000, ploidy = 2, population = 1700000, mutationrate = 0.0000000055, percentmissing = 0, percentsitemissing = 0, outputfile = 'myvcftest', samp_num = 20, sample_names = None):
    
    # when custom sample names are provided we set sample size from the names
    if sample_names is not None:
        try:
            names_len = len(sample_names)
        except Exception as e:
            print('Error, samples must be a list of names')
            sys.exit(0)
        if names_len < 1:
            print('Error, samples list must have at least one name')
            sys.exit(0)
        samp_num = names_len

    for x in range(amountofruns):
        if outputfile is not None:
            outputfilename = outputfile + str(seed) + '.vcf'
        else:
            outputfilename = 'None'
            
        sim = MyVcfSim(chrom, sitesize, ploidy, population, mutationrate, percentmissing, percentsitemissing, seed, outputfilename, 
                       samp_num, 'population.txt', 'vcf', sample_names)

        sim.simulate_vcfs()
        seed+=1

def main():
    #argument parser using argparse
    
    parser = argparse.ArgumentParser(description = 'VCFSim is used to accurately simulate VCFs')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('--chromosome', type=str, nargs = '?', help = 'Chromosome name', required = False)
    optional.add_argument('--replicates', type=int, nargs = '?', help = 'Amount of times for Simulator to run', required = False)
    required.add_argument('--seed', type=int, nargs = '?', help = 'Random seed for VCFSim to use', required = True)
    optional.add_argument('--sequence_length', type=int, nargs = '?', help = 'Size of your site', required = False)
    optional.add_argument('--ploidy', type=int, nargs = '?', help = 'Ploidy for your VCF', required = False)
    optional.add_argument('--Ne', type=int, nargs = '?', help = 'Effective population size of the simulated population', required = False)
    optional.add_argument('--mu',type=float, nargs = '?', help = 'Mutation rate in the simulated population', required = False)  
    required.add_argument('--percent_missing_sites', type=int, nargs = '?', help = 'Percent of rows missing from your VCF', required = True) #percent missing sites   
    required.add_argument('--percent_missing_genotypes', type=int, nargs = '?', help = 'Percent of samples missing from your VCF', required = True) #percent missing genotypes, percent of genotypes that are missing from samples    
    optional.add_argument('--output_file', nargs = '?', help = 'Filename of outputed vcf, will automatically be followed by seed', required = False)    

    optional.add_argument('--param_file', nargs = '?', help = 'Specified file for multiple chromosome inputs', required = False)

    #choose one way to specify samples
    #either a numeric sample size or explicit names or a file of names
    sample_group = parser.add_mutually_exclusive_group(required = True)
    sample_group.add_argument('--sample_size', type = int, nargs = '?', help = 'Amount of samples from population in VCF')
    sample_group.add_argument('--samples', nargs = '*', help = 'Custom sample names space separated')
    sample_group.add_argument('--samples_file', nargs = '?', help = 'File with one whitespace separated line of sample names')

    args = parser.parse_args()

    # read custom names if provided
    custom_names = None
    if args.samples_file is not None:
        try:
            with open(args.samples_file, 'r') as sf:
                line = sf.readline()
            tokens = line.split()
            if len(tokens) < 1:
                print('Error, samples file must contain at least one name')
                sys.exit(0)
            custom_names = tokens
        except Exception as e:
            print('Error, could not read samples file')
            sys.exit(0)
    elif args.samples is not None and len(args.samples) > 0:
        custom_names = args.samples

    if args.param_file is not None:
        # when using multiple chromosome file we forward the custom names too
        # sample size must be determined here so we never pass None
        if custom_names is not None:
            effective_samp_num = len(custom_names)
        else:
            effective_samp_num = args.sample_size

        multiple_chrom(chromfilename = args.param_file, seed = args.seed, percentmissing = args.percent_missing_sites, percentsitemissing = args.percent_missing_genotypes, outputfile = args.output_file, samp_num = effective_samp_num, sample_names = custom_names)
    
    elif (args.param_file is None and (args.chromosome is None or args.replicates is None or args.sequence_length is None
                                      or args.ploidy is None or args.Ne is None or args.mu is None)):
        print("Error, no parameter file is specified, or missing one of the following arguments: chromosome, replicates, sequence_length, ploidy, population size, or mutation rate")

    elif (args.param_file is not None and (args.chromosome is not None or args.replicates is not None or 
                                           args.sequence_length is not None or args.ploidy is not 
                                           None or args.Ne is not None or args.mu is not None)):
        print("Error,  parameter file is specified, however so are one of the following: chromosome, replicates, sequence_length, ploidy, population size, or mutation rate")
    
    elif args.mu >= 1:
        print("Error: Mutation rate must be less than 1")

    elif type(args.output_file) != str and args.output_file is not None:
        print("Error: Output_file must be a string")
        
    else:
        vcf_simulator(chrom = args.chromosome , amountofruns = args.replicates, seed = args.seed, sitesize = args.sequence_length, 
                     ploidy = args.ploidy, population = args.Ne, mutationrate = args.mu, percentmissing = args.percent_missing_sites,
                     percentsitemissing = args.percent_missing_genotypes, outputfile = args.output_file, samp_num = args.sample_size, sample_names = custom_names)
    
    #argument checker after to double check
    
warnings.filterwarnings("ignore")
main()