import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
import argparse
from IPython.display import SVG, display
from src.SimulatorClass import MyVcfSim


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
        sim.make_population_file()
        sim.simulate_vcfs()
    
        command1 = 'bgzip ' + sim.outputfile
        command2 = 'tabix -f -p vcf -0 ' + sim.outputfile + '.gz'
        command3 = 'pixy --stats pi --vcf ' + sim.outputfile + '.gz --populations '+ sim.samp_file + ' --window_size 10000 --output_prefix PixyTest' + str(sim.randoseed) 
        command4 = sim.outputfile + '.gz'
        command5 = sim.outputfile + '.gz.tbi'
        command6 = sim.samp_file
        
        #Commands used to make text files in the proper format to be inputed into pixy. Deletes excess files aswell. 
        
        os.system(command1)
        os.system(command2)
        os.system(command3)
        os.remove(command4)
        os.remove(command5)
        os.remove(command6)
        #Runs simulation with the data from above.
        
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

mass_simulator(10000, 2, 1700000, 0.0000000055, 0, 0, 1234, 'my.vcf', 20, 'population.txt', 'PixyFolderTest', 10)
#Mass simulator function to run everything above.s
#Parameters in order: Site Size, Ploidy, Population, Mutation Rate, Percent Of Missing data to be taken out,-
#Starting seed, Output VCF file, Sample Size, Sample file, Output Folder, Amount of times to run everything

#seperate pixy stuff
