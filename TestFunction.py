import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
import argparse
from IPython.display import SVG, display
from SimulatorClass import MyVcfSim

def mass_simulator(sitesize, ploidy, popnum, mutationrate,
                   percentmissing, percentsitemissing, seed, outputfile, sampnum, sampfile, foldername, amountofruns = 1):
    
    columns = ('pop_num', 'ploidy', 'mutation_rate', 'percent_missing', 'percent_site_missing', 'seed')
    df = pd.DataFrame(index = columns).T #Creates a pandas data frame with the columns above.
    
    defaultdirectory = os.getcwd() 
    os.mkdir(foldername) 
    os.chdir(foldername) 

    for x in range(amountofruns):
        
        foldername = 'folder' + str(seed)
        sim = MyVcfSim(sitesize, ploidy, popnum, mutationrate, percentmissing, percentsitemissing, seed, 
                       outputfile, sampnum, sampfile, foldername)
        
        #Runs simulation with the data from above.
        
        sim.mass_simulate_pix()
                
        pixdata = pd.read_csv('PixyTest' + str(seed) + '_pi.txt', sep="\t")
        
        #Uses pandas to read data given from pixy into another data frame
        
        data = {
        0: {'pop_num': sim.pop_num, 'ploidy': sim.ploidy, 'mutation_rate': sim.mutationrate, 'percent_missing': sim.percentmissing, 
                          'percent_site_missing': sim.percentsitemissing, 'seed': seed},
        }
        
        df2 = pd.DataFrame(data=data, index = columns).T
        
        df2 = pd.concat([df2,pixdata], axis = 1) 
        
        df = pd.concat([df, df2], axis = 0)
        
        #Combines the two data frames together side by side so all the pixy data matches up to the correct test run in a big data frame
        
        seed+=1
        
    os.chdir(defaultdirectory) 
    df.to_csv('Results.csv', index = False)
    
    #Outputs all the data from the main frame into its own csv file
    
    
#Quick function to test if VCF files are getting outputed correctly
def vcf_simulator(amountofruns = 1, seed = 1234, foldername = 'PixyFolder', sitesize = 10000, ploidy = 2, population = 1700000, 
                 mutationrate = 0.0000000055, percentmissing = 0, percentsitemissing = 0, outputfile = 'myvcftest', samp_num = 20):
    
    for x in range(amountofruns):
        outputfilename = outputfile + str(seed) + '.txt'
        sim = MyVcfSim(sitesize, ploidy, population, mutationrate, percentmissing, percentsitemissing, seed, outputfilename, samp_num, 'population.txt', 'vcf')
        sim.simulate_vcfs()
        seed+=1

#mass_simulator(10000, 2, 1700000, 0.0000000055, 0, 1234, 'my.vcf', 20, 'population.txt', 'PixyFolder', 100)
#Mass simulator function to run everything above.s
#Parameters in order: Site Size, Ploidy, Population, Mutation Rate, Percent Of Missing data to be taken out,-
#Starting seed, Output VCF file, Sample Size, Sample file, Output Folder, Amount of times to run everything

#seperate pixy stuff

def main():
    #argument parser using argparse
    
    parser = argparse.ArgumentParser(description = 'VCFSim is used to accurately simulate VCFs')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    
    required.add_argument('--replicates', type=int, nargs = '?', help = 'Amount of times for Simulator to run', required = True)
    required.add_argument('--seed', type=int, nargs = '?', help = 'Random seed for VCFSim to use', required = True)
    required.add_argument('--sequence_length', type=int, nargs = '?', help = 'Size of your site', required = True)
    required.add_argument('--ploidy', type=int, nargs = '?', help = 'Ploidy for your VCF', required = True)
    required.add_argument('--Ne', type=int, nargs = '?', help = 'Effective population size of the simulated population', required = True)
    required.add_argument('--mu',type=float, nargs = '?', help = 'Mutation rate in the simulated population', required = True)  
    required.add_argument('--percent_missing_sites', type=int, nargs = '?', help = 'Percent of rows missing from your VCF', required = True) #percent missing sites   
    required.add_argument('--percent_missing_genotypes', type=int, nargs = '?', help = 'Percent of samples missing from your VCF', required = True) #percent missing genotypes, percent of genotypes that are missing from samples    
    required.add_argument('--output_file', nargs = '?', help = 'Filename of outputed vcf, will automatically be followed by seed', required = True)    
    required.add_argument('--sample_size', type=int, nargs = '?', help = 'Amount of samples from population in VCF', required = True)    

    args = parser.parse_args()


    vcf_simulator(amountofruns = args.replicates, seed = args.seed, sitesize = args.sequence_length, ploidy = args.ploidy, 
                 population = args.Ne, mutationrate = args.mu, percentmissing = args.percent_missing_sites,
                 percentsitemissing = args.percent_missing_genotypes, outputfile = args.output_file, samp_num = args.sample_size)
    
    #argument checker after to double check

main()