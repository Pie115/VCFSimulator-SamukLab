import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
from IPython.display import SVG, display
import io

class MyVcfSim:
    
    def __init__(self, chrom, site_size, ploidy, pop_num, mutationrate, percentmissing, percentsitemissing, randoseed, outputfile, samp_num, samp_file, folder, sample_names = None):

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

        # store custom names if provided
        if sample_names is not None:
            self.sample_names = list(sample_names)
        else:
            self.sample_names = None

        # prebuild small strings used many times
        change_list_zero = []
        for i in range(self.ploidy):
            change_list_zero.append('0')
        joiner_temp = '|'
        self.change_zero = joiner_temp.join(change_list_zero)

        change_list_missing = []
        for i in range(self.ploidy):
            change_list_missing.append('.')
        self.change_missing = joiner_temp.join(change_list_missing)

        self.col_start = None
        self.col_end = None
    
    def make_site_mask(self):
        temp_percent = self.percentmissing / 100
        temp_var = self.site_size * temp_percent

        temp_array = np.zeros(self.site_size, dtype=int)

        k = round(temp_var)
        if k > 0:
            order = np.random.permutation(self.site_size)
            missing_idx = order[:k]
            temp_array[missing_idx] = 1

        return temp_array
    
    def row_changes(self, row, vcfdata, tempvcf):
        uniquelist = []

        col_start = self.col_start
        col_end = self.col_end
        altlist = row[col_start:col_end+1].values
        
        randomsitemissing = round((self.percentsitemissing / 100) * (self.samp_num - 1))

        randomsites = np.arange(1, self.samp_num)
        np.random.shuffle(randomsites)
        randomsites = randomsites[:randomsitemissing]

        change = self.change_zero
        
        for idx_site in randomsites:
            altlist[idx_site] = change
        
        refindex = 0

        if ((1 in randomsites) and ((self.percentsitemissing / 100) != 1)):
            for i in range(1, self.samp_num):
                if i not in randomsites:
                    refindex = i * self.ploidy
                    break
        
        if (self.ploidy != 1):
            for item in altlist:
                new_items = item.split('|')
                for x in new_items:
                    uniquelist.append(int(x))
        elif (self.ploidy == 1):
            for x in altlist:
                uniquelist.append(int(x))
            
        a = 'ALT'
        alt_chars = row[a]
        alt_chars = alt_chars.replace(',', '')
        alt_chars = list(alt_chars)

        referencegenome = uniquelist[refindex]
        if (referencegenome != 0):
            oldref = row['REF']
            row['REF'] = alt_chars[referencegenome - 1]
            removed_alt = alt_chars[referencegenome - 1]
            alt_chars.remove(removed_alt)
            alt_chars.append(oldref)

            for i in range(len(uniquelist)):
                if (uniquelist[i] == 0):
                    uniquelist[i] = len(alt_chars)
                elif (uniquelist[i] == referencegenome):
                    uniquelist[i] = 0
                elif (uniquelist[i] != 1):
                    uniquelist[i] = uniquelist[i] - 1
        
        uniquelist = np.array(uniquelist, dtype=int)
                    
        finaloutput = uniquelist
        unique_vals = np.unique(uniquelist)
        
        templist = []
        for i in range(len(alt_chars)):
            idx_val = i + 1
            if (idx_val in unique_vals):
                templist.append(alt_chars[i])
        alt_chars = templist

        if (len(alt_chars) == 0):
            alt_str = "."
        else:
            alt_str = ','.join(alt_chars)

        row[a] = alt_str

        if ((1 not in finaloutput) and (finaloutput.sum != 0)):
            for i in range(len(finaloutput)):
                if (finaloutput[i] > 1):
                    finaloutput[i] = finaloutput[i] - 1
                
        finaldata = []
        for i in range(0, len(finaloutput), self.ploidy):
            finaldata.append(finaloutput[i:i + self.ploidy])
        finaldataoutput = []
        for data in finaldata:
            finaldataoutput.append("|".join(str(i) for i in data))
        
        change_missing_local = self.change_missing
        for idx_site in randomsites:
            finaldataoutput[idx_site] = change_missing_local
        
        if ((self.percentsitemissing / 100) == 1):
            row['REF'] = '.'
        
        row[col_start:col_end+1] = finaldataoutput

        return row


    def make_missing_vcf(self, ts):
        np.random.seed(self.randoseed)
        site_mask = self.make_site_mask()

        buf = io.StringIO()
        ts.write_vcf(buf, site_mask=site_mask)
        buf.seek(0)

        header_lines = []
        data_lines = []

        linenum = 1
        for line in buf:
            if (linenum <= 5):
                header_lines.append(line)
            else:
                if (linenum == 6):
                    line = line[1:]
                data_lines.append(line)
            linenum = linenum + 1

        body_buf = io.StringIO()
        for line in data_lines:
            body_buf.write(line)
        body_buf.seek(0)

        vcfdata = pd.read_csv(body_buf, delimiter = '\t')
        
        a = 'tsk_0'
        self.col_start = vcfdata.columns.get_loc(a)
        self.col_end = vcfdata.columns.get_loc(f"tsk_{self.samp_num - 1}")
        
        tempvcf = vcfdata
        
        vcfdata = vcfdata.apply(self.row_changes, axis = 1, args = (vcfdata, tempvcf))
 
        vcfdata["CHROM"] = self.chrom
        vcfdata["POS"] = vcfdata["POS"] + 1
        vcfdata["ID"] = '.'

        if ("tsk_0" in vcfdata.columns):
            del vcfdata["tsk_0"]

        # rename tsk style columns to custom names when provided
        if (self.sample_names is not None):
            # we rename tsk_1 to first name and so on
            idx_val = 1
            name_index = 0
            while(idx_val < self.samp_num and name_index < len(self.sample_names)):
                old_name = "tsk_" + str(idx_val)
                new_name = self.sample_names[name_index]
                if old_name in vcfdata.columns:
                    vcfdata.rename(columns={old_name: new_name}, inplace=True)
                idx_val = idx_val + 1
                name_index = name_index + 1
        
        if (self.outputfile != 'None'):
            with open(self.outputfile, 'w') as fout:

                for line in header_lines:
                    fout.write(line)

                #fout.write('##vcfsim_version=1.0\n')
                #fout.write('##vcfsim_command=' + ' '.join(sys.argv) + '\n')

                csv_buf = io.StringIO()
                vcfdata.to_csv(csv_buf, index=False, sep='\t', header=True)
                csv_text = csv_buf.getvalue()
                fout.write('#')
                fout.write(csv_text)
        else:
            display(vcfdata.to_string())
        
    def make_population_file(self):
        np.random.seed(self.randoseed)
        file = open(self.samp_file, "a")

        for x in range(0, self.samp_num):
            file.write("tsk_"+ str(x) + "\t1\n")
        file.close()

    def simulate_vcfs(self):
        np.random.seed(self.randoseed)
        ts = msprime.sim_ancestry(samples=[msprime.SampleSet(self.samp_num, ploidy=self.ploidy)], population_size=self.pop_num, random_seed=self.randoseed, sequence_length=self.site_size)

        rlg = np.random.default_rng(self.randoseed)
        rchar = rlg.integers(low=0, high=4)

        difference_counter = range(self.site_size)

        tables = ts.dump_tables()

        for x in difference_counter:
            tables.sites.add_row(x, "A")

        ts = tables.tree_sequence()
        ts = msprime.sim_mutations(ts, rate=self.mutationrate, random_seed=self.randoseed)

        self.make_missing_vcf(ts)
