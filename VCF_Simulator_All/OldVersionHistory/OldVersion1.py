import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
from IPython.display import SVG, display

class MyVcfSim:
    
    def __init__(self, site_size, ploidy, pop_num, mutationrate, percentmissing, randoseed, outputfile, samp_num, samp_file, folder):
        
        self.site_size = site_size
        self.ploidy = ploidy
        self.pop_num = pop_num
        self.mutationrate = mutationrate
        self.percentmissing = percentmissing
        self.randoseed = randoseed
        self.outputfile = outputfile
        self.samp_num = samp_num
        self.samp_file = samp_file
        self.folder = folder
    
    def make_site_mask(self):

        temp_percent = self.percentmissing/100 #takes percent of user input and converts it to actual percentage
        temp_var = self.site_size*temp_percent  #finds amount of VCF's that need to be taken out
        counter = range(int(temp_var)) #Makes a counter to traverse VCF's
        temp_array = np.zeros(self.site_size, dtype = int) #Makes array of 0's 
        
        mylist = np.arange(0, self.site_size,1) #Generates an array of 1-size_size
        rand_array = np.random.choice(mylist, self.site_size, replace=False) #fills a new array with random non repeating values from last step
        

        print('Total Amount:', temp_array.size)  #Shows total amount of sites

        print('----------------')

        for x in counter:
            temp_array[rand_array[x]] = 1 #Changes value to 1 which deletes the value

        print('Amount Missing:', np.sum(temp_array))

        print('----------------')
        
        rand_array = np.empty(self.site_size, dtype = int)
                
        return temp_array

    def make_missing_vcf(self, ts):

        site_mask = self.make_site_mask() #Makes sites from the number of sites in ts

        with open(self.outputfile, "w") as f:
            ts.write_vcf(f, site_mask=site_mask)

    def make_population_file(self):

        file = open(self.samp_file, "a") #Makes a text file full of population data as a parameter for pixy

        for x in range(0, self.samp_num):
            file.write("tsk_"+ str(x) + "\t1\n")
        file.close()

    def simulate_vcfs(self):
        
        ts = msprime.sim_ancestry(samples=[msprime.SampleSet(self.samp_num, ploidy=self.ploidy)], population_size = self.pop_num, random_seed=self.randoseed, sequence_length = self.site_size)

        #sequence length and site size have to be set to SAME number

        #sets number of sites

        rlg = np.random.default_rng(self.randoseed) #random letter generator 

        rchar = rlg.integers(low=0, high=4) #randomly choses number from 0 to 3 which represents genome

        difference_counter = range(self.site_size) 

        tables = ts.dump_tables() 

        for x in difference_counter: #Randomly assings reference genome for each site
            
            tables.sites.add_row(x, "A") 

        ts = tables.tree_sequence()

        ts = msprime.sim_mutations(ts, rate=self.mutationrate, random_seed=self.randoseed) #Mutates the sites at random

        self.make_missing_vcf(ts)

        #os.system('cat my.vcf | grep "^#" > vcf_new.vcf')
        #os.system('cat my.vcf | grep -v "^#" | awk -v s=1 \'{$2=$2+s; $3=$3+s; print}\' >> vcf_new.vcf')
        #os.system('rm my.vcf')

    def simulate_pix_once(self):
        #Function meant to test pixy once, keeping all the files
        defaultdirectory = os.getcwd() 

        newfolder = self.folder

        os.mkdir(newfolder) 
        os.chdir(newfolder) #Creates new directory to store data files

        self.make_population_file()
        self.simulate_vcfs() 

        command1 = 'bgzip ' + self.outputfile
        command2 = 'tabix -f -p vcf -0 ' + self.outputfile + '.gz'
        command3 = 'pixy --stats pi --vcf ' + self.outputfile + '.gz --populations '+ self.samp_file + ' --window_size 10000 --output_prefix PixyTest' + str(self.randoseed) 

        #Commands used to make text files in the proper format to be inputed into pixy
        
        os.system(command1)
        os.system(command2)
        os.system(command3)
        
        os.chdir(defaultdirectory) 
            
    def mass_simulate_pix(self):
        #Function meant to test pixy multiple times. All files except for the main pixy file are deleted. 
        self.make_population_file()
        self.simulate_vcfs()

        command1 = 'bgzip ' + self.outputfile
        command2 = 'tabix -f -p vcf -0 ' + self.outputfile + '.gz'
        command3 = 'pixy --stats pi --vcf ' + self.outputfile + '.gz --populations '+ self.samp_file + ' --window_size 10000 --output_prefix PixyTest' + str(self.randoseed) 
        command4 = self.outputfile + '.gz'
        command5 = self.outputfile + '.gz.tbi'
        command6 = self.samp_file
        
        #Commands used to make text files in the proper format to be inputed into pixy. Deletes excess files aswell. 
        
        os.system(command1)
        os.system(command2)
        os.system(command3)
        os.remove(command4)
        os.remove(command5)
        os.remove(command6)
        
#Function meant to simulate pixy on a mass scale, outputting all data to a .csv data sheet
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
    df.to_csv('data.csv')
    
    #Outputs all the data from the main frame into its own csv file
    
    
#Quick function to test if VCF files are getting outputed correctly
def vcf_simulator(amountofruns = 1, seed = 1000, foldername = 'PixyFolder'):
    
    for x in range(amountofruns):
        outputfile = 'myvcf' + str(seed)
        foldername = 'folder' + str(seed)
        sim = MyVcfSim(10, 1, 1700000, 0.0000001, 20, seed, outputfile, 20, 'population.txt', 'vcff')
        sim.simulate_vcfs()
        seed+=1

mass_simulator(10000, 2, 1700000, 0.0000001, 20, 1000, 'my.vcf', 20, 'population.txt', 'PixyFolder', 10)
#Mass simulator function to run everything above.
#Parameters in order: Site Size, Ploidy, Population, Mutation Rate, Percent Of Missing data to be taken out,-
#Starting seed, Output VCF file, Sample Size, Sample file, Output Folder, Amount of times to run everything