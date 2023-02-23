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
            
        f = open(self.outputfile, 'r')
        d = open('mytempvcf.txt', 'w')

        a = f.readline()
        linenum = 1

        while(a):
            if(linenum > 5):
                if(linenum == 6):
                    a = a[1:]
                d.write(a)
            a = f.readline()
            linenum += 1

        f.close()
        d.close()

        vcfdata = pd.read_csv('mytempvcf.txt', delimiter = '\t')
        tempvcf = pd.read_csv('mytempvcf.txt', delimiter = '\t')

        a = 'tsk_0'
        
        for i in range(self.samp_num):
            a = a.replace(str(i-1), str(i))
            if(self.ploidy != 1):
                tempvcf[a] = tempvcf[a].str.replace(r'|', '')

        rows = len(tempvcf.index)
        iteration = 0

        while(iteration < rows):
            uniquelist = []
            a = 'tsk_0'

            for i in range(self.samp_num):
                a = a.replace(str(i-1), str(i))
                altlist = tempvcf[a][iteration]
                #print(altlist)
                if(self.ploidy != 1):
                    uniquelist.extend(eval(i) for i in list(altlist))
                elif(self.ploidy == 1):
                    uniquelist.append(tempvcf[a][iteration])

            a = 'ALT'
            altlist = tempvcf[a][iteration]
            altlist = altlist.replace(r',','')
            altlist = list(altlist)
            
            #print(altlist)
            uniquelist = np.unique(uniquelist)
            #print(uniquelist)
            if(len(uniquelist) == 1):
                while(len(altlist) > len(uniquelist)):
                    altlist.pop()
            elif(sum(uniquelist) == 0):
                while(len(altlist) > 0):
                    altlist.pop()
            else:
                while(len(altlist) >= len(uniquelist)):
                    altlist.pop()

            if(len(altlist) == 0):
                altlist = "."
            else:
                altlist = ','.join(altlist)

            vcfdata[a][iteration] = altlist
            #print(altlist)
            #print(uniquelist)
            #print("\n")
            iteration += 1

        #print(tempvcf)
        #print(vcfdata)

        f = open(self.outputfile, 'r')
        linenum = 1
        filearr = []
        a = f.readline()

        while(a):
            if(linenum < 6):
                filearr.append(a)
            a = f.readline()
            linenum += 1
        os.remove(self.outputfile)

        vcfdata.to_csv(self.outputfile, mode = 'a', index = False, sep = '\t', header = True)

        with open(self.outputfile, 'r+') as f: 
            tempdata = f.read() 
            f.seek(0, 0) 
            for lines in filearr:
                f.write(lines)
            f.write('#' + tempdata)

        #os.remove('mytempvcf.txt')
        
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
        