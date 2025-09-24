<img align="right" width="160" src="https://github.com/user-attachments/assets/228cfba0-bec0-4b74-8010-412d0f184417">

# vcfsim
vcfsim is a new command-line tool for generating simulated VCF's (variant call format files for encoding genetic data). Leveraging a coalescent simulating backend and providing an interface from Msprime coalescent simulating package to pandas. VCF's can now be easily simulated with just a few command line arguments!

## Author 
Paimon Goulart (UC Riverside)
Kieran Samuk (UC Riverside)

## Installation
First create and activate a conda environment for vcfsim:

```shell
conda create -n vcfsim_env python=3.10
conda activate vcfsim_env
```

vcfsim is currently available on bioconda, and can be installed by using the following command:
```shell
conda install bioconda::vcfsim
```

For more detailed installation instructions, please visit:  
https://bioconda.github.io/recipes/vcfsim/README.html?highlight=vcfsi#package-package%20&#x27;vcfsim&#x27;

## Arguments 
Here is the list of required/optional arguments to run vcfsim

### Required
--seed [SEED] Random seed for vcfsim to use  

--percent_missing_sites [PERCENT_MISSING_SITES] Percent of rows missing from your VCF  

--percent_missing_genotypes [PERCENT_MISSING_GENOTYPES] Percent of samples missing from your VCF  

One of the following three options must also be provided to set the samples:  
- --sample_size [SAMPLE_SIZE] Amount of samples from population in VCF  
- --samples [SAMPLES ...] Custom sample names, space separated (e.g. A1 B1 C1)  
- --samples_file [SAMPLES_FILE] File containing one whitespace separated line of custom sample names  

### Optional
--chromosome [CHROMOSOME] Chromosome name  

--replicates [REPLICATES] Amount of times for Simulator to run  

--sequence_length [SEQUENCE_LENGTH] Size of your site  

--ploidy [PLOIDY] Ploidy for your VCF  

--Ne [NE] Effective population size of the simulated population  

--mu [MU] Mutation rate in the simulated population  

--output_file [OUTPUT_FILE] Filename of outputed vcf, will automatically be followed by seed  

--chromosome_file [CHROMOSOME_FILE] Specified file for multiple chromosome inputs  

--population_mode [1|2] Mode of population structure to simulate. 1 = single population (default), 2 = structured population (C splits into A & B).

--time [TIME] Split time for population mode 2 (e.g. generations before present). Required if --population_mode 2 is specified.

## Usage
Typical usage for vcfsim can be done by using the following command:  

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --sample_size 10
```

This will create a vcf output file by the name of myvcf1234, or myvcf followed by the seed given for the input.  
If input for replicates was given as a higher number than 1, 2 for example, then vcfsim will create two output files by the name of myvcf1234 and myvcf1235, adding one to the seed after every run.  

NOTE: An output file doesn't needed to be specified. If no output file is specified, then the vcf will be outputed to the command line.

Screenshot of output file:
<img width="1437" height="458" alt="Image" src="https://github.com/user-attachments/assets/11078b68-6a62-44e0-bf8d-87c34544b2a6" />

### Using custom sample names
Instead of `--sample_size`, you can provide explicit names:  

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --samples A1 B1 C1 D1
```

This will automatically set the sample size to 4 and label the VCF columns `A1 B1 C1 D1`.  

You can also read the names from a file containing a single whitespace separated line:  

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --samples_file names.txt
```

Where `names.txt` might contain:  
```
A1 B1 C1 D1 E1
```

Otherwise, sample identifiers will default to tsk_0,...,tsk_n


### Simulating a structured population split (population_mode = 2)
To simulate a demographic split between populations A and B from an ancestral population C:

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --sample_size 10 --population_mode 2 --time 1000
```

### Multiple chromosome inputs
Another way vcfsim can be used is by providing a file for multiple chromosome inputs.  

Your input file should be in the form of a text file, and should be formatted as such:  
<img width="221" alt="Example input file" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/39ca4b31-c58e-4fff-8b52-456849678339">

The columns are in the order of: chromosone, ploidy, sequence length, population size, mutation rate.  
Each row will represent a seperate run of vcfsim, all these runs will be concatenated to the same file in the end.  

The following command should be used when running vcfsim in this way:

```shell
vcfsim --seed 1234 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --sample_size 10 --chromosome_file input.txt
```

You can also combine a param file with custom names:

```shell
vcfsim --seed 1234 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --samples_file names.txt --chromosome_file input.txt
```

When done this way, the output should look like such:  
<img width="1437" height="458" alt="Image" src="https://github.com/user-attachments/assets/11078b68-6a62-44e0-bf8d-87c34544b2a6" />

With the concatenated vcf looking like:  
<img width="1059" alt="ExampleInput" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/fb6508eb-34cd-473a-bc19-762858ed4c31">
