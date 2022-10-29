import msprime
import numpy as np
import sys
import os
import time
from IPython.display import SVG, display

def make_site_mask(size, percentmissing = 0, seed=0):
    
    temp_percent = percentmissing/100 #takes percent of user input and converts it to actual percentage
    temp_var = size*temp_percent  #finds amount of VCF's that need to be taken out
    counter = range(int(temp_var)) #Makes a counter to traverse VCF's
    temp_array = np.zeros(size, dtype = int) #Makes array of 0's 
    rng = np.random.default_rng(seed) #Makes random number
    rints = rng.integers(low=0, high=size) #Generates random number between 0 and the amount of sites
    
    rand_array = np.empty(size, dtype = int) #Makes empty array to store random variables
    
    
    print('Total Amount:', temp_array.size)  #Shows total amount of sites
    
    print('----------------')
    
    for x in counter:
        rints = rng.integers(low=0, high=size) #Generates random number between 0 and the amount of sites
        while rints in rand_array:
            rints = rng.integers(low=0, high=size) #Looks through array to find if the random generated number has already been repeated
        rand_array[x] = rints #stores random number in array
        temp_array[rints] = 1 #Changes value to 1 which deletes the value

    print('Amount Missing:', np.sum(temp_array))
    
    print('----------------')
    
    return temp_array

def make_missing_vcf(ts, percentmissing = 0, outputfile = ''):
    
    site_mask = make_site_mask(ts.num_sites, percentmissing, 1234) #Makes sites from the number of sites in ts
    
    with open(outputfile, "w") as f:
        ts.write_vcf(f, site_mask=site_mask)

def make_population_file(pop_num = 20, pop_file = ''):

    file = open(pop_file, "a")

    for x in range(0, pop_num):
        file.write("tsk_"+ str(x) + "\t1\n")
    file.close()
    
def simulate_vcfs(site_size = 10000, percentmissing = 0, outputfile = '', pop_num = 20):
    ts = msprime.sim_ancestry(pop_num, random_seed=1, sequence_length = site_size)

    #sequence length and site size have to be set to SAME number

    #sets number of sites
    
    rlg = np.random.default_rng(1234) #random letter generator 

    rchar = rlg.integers(low=0, high=4) #randomly choses number from 0 to 3 which represents genome

    difference_counter = range(site_size) 

    tables = ts.dump_tables() 

    for x in difference_counter: #Randomly assings reference genome for each site
    
        rchar = rlg.integers(low=0, high=4)
        
        if rchar == 0:
            letter = "A"
        
        elif rchar == 1:
            letter = "C"
        
        elif rchar == 2:
            letter = "G"
        
        elif rchar == 3:
            letter = "T"
        
        tables.sites.add_row(x, "A") 
    
    ts = tables.tree_sequence()

    ts = msprime.sim_mutations(ts, rate=0.1, random_seed=1234) #Mutates the sites at random

    make_missing_vcf(ts, percentmissing, outputfile)

    #os.system('cat my.vcf | grep "^#" > vcf_new.vcf')
    #os.system('cat my.vcf | grep -v "^#" | awk -v s=1 \'{$2=$2+s; $3=$3+s; print}\' >> vcf_new.vcf')
    #os.system('rm my.vcf')
    
def simulate_pix(site_size = 10000, percentmissing = 0, outputfile = '', pop_num = 20, pop_file = '', folder = ''):
    
    defaultdirectory = os.getcwd()
    
    newfolder = folder
    
    os.mkdir(newfolder) 
    os.chdir(newfolder) 
    
    make_population_file(pop_num, pop_file)
    simulate_vcfs(site_size, percentmissing, outputfile, pop_num)
    
    command1 = 'bgzip ' + outputfile
    command2 = 'tabix -f -p vcf -0 ' + outputfile + '.gz'
    command3 = 'pixy --stats pi --vcf ' + outputfile + '.gz --populations '+ pop_file + ' --window_size 10000'
    
    
    os.system(command1)
    os.system(command2)
    os.system(command3)
    
    os.chdir(defaultdirectory) 
    
simulate_pix(10000, 20, 'my.vcf', 20, 'population.txt', 'myvcfholder')