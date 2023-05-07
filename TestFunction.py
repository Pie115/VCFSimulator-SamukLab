import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
from IPython.display import SVG, display
from SimulatorClass import MyVcfSim

def mass_simulator(sitesize, ploidy, popnum, mutationrate,
                   percentmissing, seed, outputfile, sampnum, sampfile, foldername, amountofruns = 1):
    
    columns = ('ploidy', 'mutation_rate', 'percent_missing', 'seed')
    df = pd.DataFrame(index = columns).T #Creates a pandas data frame with the columns above.
    
    defaultdirectory = os.getcwd() 
    os.mkdir(foldername) 
    os.chdir(foldername) 

    for x in range(amountofruns):
        
        foldername = 'folder' + str(seed)
        sim = MyVcfSim(sitesize, ploidy, popnum, mutationrate, percentmissing, seed, 
                       outputfile, sampnum, sampfile, foldername)
        
        #Runs simulation with the data from above.
        
        sim.mass_simulate_pix()
                
        pixdata = pd.read_csv('PixyTest' + str(seed) + '_pi.txt', sep="\t")
        
        #Uses pandas to read data given from pixy into another data frame
        
        data = {
        0: {'ploidy': sim.ploidy, 'mutation_rate': sim.mutationrate, 'percent_missing': sim.percentmissing, 
                          'seed': seed},
        }
        
        df2 = pd.DataFrame(data=data, index = columns).T
        
        df2 = pd.concat([df2,pixdata], axis = 1) 
        
        df = pd.concat([df, df2], axis = 0)
        
        #Combines the two data frames together side by side so all the pixy data matches up to the correct test run in a big data frame
        
        seed+=1
        
    os.chdir(defaultdirectory) 
    df.to_csv('CheckerResults.csv', index = False)
    
    #Outputs all the data from the main frame into its own csv file
    
    
#Quick function to test if VCF files are getting outputed correctly
def vcf_simulator(amountofruns = 1, seed = 1234, foldername = 'PixyFolder'):
    
    for x in range(amountofruns):
        outputfile = 'myvcf' + str(seed) + '.txt'
        foldername = 'folder' + str(seed)
        sim = MyVcfSim(10000, 3, 1700000, 0.0000000055, 0, seed, outputfile, 20, 'population.txt', 'vcf')
        sim.simulate_vcfs()
        seed+=1

#mass_simulator(10000, 2, 1700000, 0.0000000055, 0, 1234, 'my.vcf', 20, 'population.txt', 'PixyFolder', 100)
#Mass simulator function to run everything above.
#Parameters in order: Site Size, Ploidy, Population, Mutation Rate, Percent Of Missing data to be taken out,-
#Starting seed, Output VCF file, Sample Size, Sample file, Output Folder, Amount of times to run everything

vcf_simulator()
