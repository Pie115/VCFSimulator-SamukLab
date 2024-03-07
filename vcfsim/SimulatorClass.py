import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
from IPython.display import SVG, display

class MyVcfSim:
    
    def __init__(self, chrom, site_size, ploidy, pop_num, mutationrate, percentmissing, percentsitemissing, randoseed, outputfile, samp_num, samp_file, folder):

        self.chrom = chrom
        self.site_size = site_size
        self.ploidy = ploidy
        self.pop_num = pop_num
        self.mutationrate = mutationrate
        self.percentmissing = percentmissing
        self.percentsitemissing = percentsitemissing
        self.randoseed = randoseed
        self.outputfile = outputfile
        self.samp_num = samp_num + 1
        self.samp_file = samp_file
        self.folder = folder
    
    def make_site_mask(self):
        #np.random.seed(self.randoseed)
        temp_percent = self.percentmissing/100 #takes percent of user input and converts it to actual percentage
        temp_var = self.site_size*temp_percent  #finds amount of VCF's that need to be taken out
        counter = range(int(temp_var)) #Makes a counter to traverse VCF's
        temp_array = np.zeros(self.site_size, dtype = int) #Makes array of 0's 
        
        mylist = np.arange(0, self.site_size,1) #Generates an array of 1-size_size
        rand_array = np.random.choice(mylist, self.site_size, replace=False) #fills a new array with random non repeating values from last step
        

        #print('Total Amount:', temp_array.size)  #Shows total amount of sites

        #print('----------------')

        for x in counter:
            temp_array[rand_array[x]] = 1 #Changes value to 1 which deletes the value

        #print('Amount Missing:', np.sum(temp_array))

        #print('----------------')
        
        rand_array = np.empty(self.site_size, dtype = int)
                
        return temp_array
    
    def row_changes(self, row, vcfdata, tempvcf):
        #np.random.seed(self.randoseed)
        uniquelist = []

        a = 'tsk_0'
        col_start = vcfdata.columns.get_loc(a)
        col_end = vcfdata.columns.get_loc(f"tsk_{self.samp_num-1}")
        altlist = row[col_start:col_end+1].values
        
        ############################################################
        randomsitemissing = (self.percentsitemissing/100) * self.samp_num
        randomsitemissing = int(randomsitemissing)
        randomsites = np.arange(0, self.samp_num)
        np.random.shuffle(randomsites)
        randomsites = randomsites[:randomsitemissing]
        #print(randomsitemissing)
        #print(randomsites)
        #print(altlist)
        change = []    #initialize default change
        for i in range(self.ploidy):
            change.append('0') 
        seperate = '|'
        change = seperate.join(change)
        
        for sites in randomsites:
            #print(altlist[sites])                         
            altlist[sites] = change
        
        refindex = 0
        
        if(0 in randomsites and (self.percentsitemissing/100) != 1):
            #print("Problem")
            for i in range(0, self.samp_num):
                if i not in randomsites:
                    refindex = i * self.ploidy
                    #print("Resolved", refindex, i)
                    break
            #print(refindex)
        #print(change)
        ############################################################
        
        if(self.ploidy != 1):
            for item in altlist:
                new_items = item.split('|')
                uniquelist.extend([int(x) for x in new_items])
        elif(self.ploidy == 1):
            uniquelist.extend([int(x) for x in altlist])
            
        a = 'ALT'
        altlist = row[a]
        altlist = altlist.replace(r',','')
        altlist = list(altlist)

        referencegenome = uniquelist[refindex]
        if(referencegenome != 0):
            oldref = row['REF']
            row['REF'] = altlist[referencegenome-1]
            altlist.remove(altlist[referencegenome-1])
            altlist.append(oldref)

            for i in range(len(uniquelist)):
                if(uniquelist[i] == 0):
                    uniquelist[i] = len(altlist)
                elif(uniquelist[i] == referencegenome):
                    uniquelist[i] = 0
                elif(uniquelist[i] != 1):
                    uniquelist[i] -= 1
        
        uniquelist = np.array(uniquelist)
                    
        finaloutput = uniquelist   #getting temp array
        uniquelist = np.unique(uniquelist)
        
        templist = []
        for i in range(len(altlist)):
            if (i+1) in uniquelist:          
                templist.append(altlist[i])
        altlist = templist

        if(len(altlist) == 0):
            altlist = "."
        else:
            altlist = ','.join(altlist)

        row[a] = altlist

        if (1 not in finaloutput and finaloutput.sum != 0): #check for edge case where array only has 2's and 0's
            for i in range(len(finaloutput)):
                if(finaloutput[i] > 1):
                    finaloutput[i] -= 1
                
        finaldata = [finaloutput[i:i+self.ploidy] for i in range(0, len(finaloutput), self.ploidy)]
        finaldataoutput = ["|".join(str(i) for i in data) for data in finaldata]
        
        ######################################
        change = []
        for i in range(self.ploidy):
            change.append('.') 
        seperate = '|'
        change = seperate.join(change)
    
        for sites in randomsites:
            #print(finaldataoutput[sites])                         
            finaldataoutput[sites] = change
        
        if(self.percentsitemissing/100 == 1):
            row['REF'] = '.'
        
        ######################################
        
        a = 'tsk_0'
        col_start = vcfdata.columns.get_loc(a)
        col_end = vcfdata.columns.get_loc(f"tsk_{self.samp_num-1}")
        row[col_start:col_end+1] = finaldataoutput

        return row


    def make_missing_vcf(self, ts):
        np.random.seed(self.randoseed)
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
        col_start = tempvcf.columns.get_loc(a)
        col_end = tempvcf.columns.get_loc(f"tsk_{self.samp_num-1}")
        cols = tempvcf.columns[col_start:col_end+1]
        tempvcf.loc[:, cols] = tempvcf.loc[:, cols].apply(lambda x: x.str.replace(r'|', '') if self.ploidy != 1 else x)
        
        rows = len(tempvcf.index)
        iteration = 0
        
        vcfdata = vcfdata.apply(self.row_changes, axis=1, args = (vcfdata, tempvcf)) #Changed while loop to a pandas apply function
        #that uses a pre compiled rowchanges function with extremely optimized 
 
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

        vcfdata["CHROM"] = self.chrom
        vcfdata["POS"] += 1
        vcfdata["ID"] = '.'

        del vcfdata["tsk_0"]
        
        if(self.outputfile != 'None'):
            vcfdata.to_csv(self.outputfile, mode = 'a', index = False, sep = '\t', header = True)
        
            with open(self.outputfile, 'r+') as f: 
                tempdata = f.read() 
                f.seek(0, 0) 
                for lines in filearr:
                    f.write(lines)
                f.write('#' + tempdata)

        else:
            display(vcfdata.to_string())

        os.remove('mytempvcf.txt')
        
    def make_population_file(self):
        np.random.seed(self.randoseed)
        file = open(self.samp_file, "a") #Makes a text file full of population data as a parameter for pixy

        for x in range(0, self.samp_num):
            file.write("tsk_"+ str(x) + "\t1\n")
        file.close()

    def simulate_vcfs(self):
        np.random.seed(self.randoseed)
        ts = msprime.sim_ancestry(samples=[msprime.SampleSet(self.samp_num, ploidy=self.ploidy)], population_size = self.pop_num, random_seed=self.randoseed, sequence_length = self.site_size)
        #sample here without changes^

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
        #with open(self.outputfile, "w") as f:
            #ts.write_vcf(f)
        self.make_missing_vcf(ts)

        #os.system('cat my.vcf | grep "^#" > vcf_new.vcf')
        #os.system('cat my.vcf | grep -v "^#" | awk -v s=1 \'{$2=$2+s; $3=$3+s; print}\' >> vcf_new.vcf')
        #os.system('rm my.vcf')