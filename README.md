# VCFSim
VCFSim is a new command-line tool for generating simulated VCF's(variant call format files for encoding genetic data). Leveraging a coalescent simulating backend and providing an interface from Msprime coalescent simulating package to pandas. VCF's can now be easily simulated with just a few command line arguments!

## Author 
Paimon Goulart (UC Riverside)

## Installation
Vcfsim is currently available on bioconda, and can be installed by using the following command:
```shell
conda install bioconda::vcfsim
```

For more detailed installation instructions, please visit:  
https://bioconda.github.io/recipes/vcfsim/README.html?highlight=vcfsi#package-package%20&#x27;vcfsim&#x27;

## Arguments 
Here is the list of required/optional arguments to run vcfsim

### Required
--seed [SEED] Random seed for VCFSim to use  

--percent_missing_sites [PERCENT_MISSING_SITES] Percent of rows missing from your VCF  

--percent_missing_genotypes [PERCENT_MISSING_GENOTYPES] Percent of samples missing from your VCF  

--sample_size [SAMPLE_SIZE] Amount of samples from population in VCF  

### Optional
--chromosome [CHROMOSOME] Chromosome name  

--replicates [REPLICATES] Amount of times for Simulator to run  

--sequence_length [SEQUENCE_LENGTH] Size of your site  

--ploidy [PLOIDY] Ploidy for your VCF  

--Ne [NE] Effective population size of the simulated population  

--mu [MU] Mutation rate in the simulated population  

--output_file [OUTPUT_FILE] Filename of outputed vcf, will automatically be followed by seed  

--param_file [PARAM_FILE] Specified file for multiple chromosome inputs  

## Usage
Typical usage for vcfsim can be done by using the following command:  

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf  --sample_size 10
```

This will create a vcf output file by the name of myvcf1234, or myvcf followed by the seed given for the input.  
If input for replicates was given as a higher number than 1, 2 for example, then vcfsim will create two output files by the name of myvcf1234 and myvcf1235, adding one to the seed after every run.  

NOTE: An output file doesn't needed to be specified. If no output file is specified, then the vcf will be outputed to the command line.

Screenshot of output file:
<img width="1065" alt="Vcfsim Screenshot" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/ec374ec4-b57f-4666-9478-cd3462c46e4f">

Another way vcfsim can be used is by providing a file for multiple chromosome inputs.  

Your input file should be in the form of a text file, and should be formatted as such:  
<img width="221" alt="Example input file" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/39ca4b31-c58e-4fff-8b52-456849678339">

The columns are in the order of: chromosone, ploidy, sequence length, population size, mutation rate.  
Each row will represent a seperate run of vcfsim, all these runs will be concatenated to the same file in the end.  

The following command should be used when running vcfsim in this way:

```shell
vcfsim  --seed 1234  --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf  --sample_size 10 --param_file input.txt
```

When done this way, the output should look like such:  
<img width="1062" alt="Inputfile" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/a16c488a-3ad0-4138-a8e0-d3abac245c7f">

With the concatenated vcf looking like:  
<img width="1059" alt="ExampleInput" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/fb6508eb-34cd-473a-bc19-762858ed4c31">

