# VCFSim
VCFSim is a new command-line tool for generating simulated VCF's(variant call format files for encoding genetic data). Leveraging a coalescent simulating backend and providing an interface from Msprime coalescent simulating package to pandas. VCF's can now be easily simulated with just a few command line arguments!

## Author 
Paimon Goulart (UC Riverside)

## Installation
Vcfsim is currently available on bioconda, and can be installed by using the following command:
```shell
conda install bioconda::vcfsim
```

For more installation instructions, please visit: 
https://bioconda.github.io/recipes/vcfsim/README.html?highlight=vcfsi#package-package%20&#x27;vcfsim&#x27;

## Arguments 
Here is the list of required/optional arguments to run vcfsim

### Required
--seed [SEED] Random seed for VCFSim to use  

--percent_missing_sites [PERCENT_MISSING_SITES] Percent of rows missing from your VCF  

--percent_missing_genotypes [PERCENT_MISSING_GENOTYPES] Percent of samples missing from your VCF  

--sample_size [SAMPLE_SIZE] Amount of samples from population in VCF  


