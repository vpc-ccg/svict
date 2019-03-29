**SViCT**: Structural Variant detrction in Circulating Tumor DNA
===================
SViCT is a computational tool for detecting structural variations from cell free DNA (cfDNA) containing low dilutions of circulating tumor DNA (ctDNA).

# Table of contents
1. [Installation](#installation)
2. [Running SViCT](#Running-SViCT)
3. [Publication](#publication)
4. [Contact & Support](#contact-support)

# Installation 

## BIOCONDA

SViCT can be istalled using [conda](https://conda.io/) package manager via [bioconda](https://bioconda.github.io/) channel:
```
To be updated
```
## Installation from Source
> *Prerequisite.* You will need g++ 4.9 and higher to compile the source code.

To install SViCT, you first need to fetch it from the [github repository](https://github.com/vpc-ccg/svict). After downloading, change the current directory to the source directory ```svict``` and run ```make``` in the terminal to create the binary file **svict**.
```
git clone https://github.com/vpc-ccg/svict.git
cd svict
make
```

If you are interested in a particular version, after downloading the git repo, checkout that version and do make.

```
git clone https://github.com/vpc-ccg/svict.git
cd svict
git checkout v1.0.0
make
```

You can also go to [releases page](https://github.com/vpc-ccg/svict/releases), click on the desired version and then click either **Source Code (zip)** or **Source code (tar.gz)** link to download the file. After decompressing it, you just switch to the ```svict``` directory and run ```make```.




# Running SViCT
SViCT requires **coordinate-sorted BAM/SAM** and **reference genome FASTA** files to run detections:

```
./svict -i input.sorted.bam -r human_genome.fa -o out
```
The output based on *input.sorted.bam* will be written to **out.vcf** in the current folder.

## Test Datasets
To grab sample data and test ```SViCT```, please first download the following two files ( ~ 3GB in total):
```
curl -L https://ndownloader.figshare.com/files/12380225 --output sim.75.sorted.bam
curl -L https://ndownloader.figshare.com/files/10144653 --output Homo_sapiens.GRCh38.87.dna.chromosomes.fa
```

Type the following command to run ```svict```:
```
./svict -i sim.75.sorted.bam -r Homo_sapiens.GRCh38.87.dna.chromosomes.fa -o test.75
```
The VCF file with the prediction results will be generated as **test.75.vcf** in the current directory.

> You need to the paths correctly for ***Homo_sapiens.GRCh38.87.dna.chromosomes.fa*** and ***sim.100.sorted.bam*** if you download these files to a folder other the ```svict``` folder. 


## Command Options ## 
### Mandatory Parameters ###
1. **-i|--input** [FILE]: Input alignment file. This file should be a sorted SAM or BAM file.
1. **-r|--reference** [FILE]: Reference geneome that the reads are algined to.

### Main Optional Parameters ###
1. **-o|--output** [STRING]: Prefix of output vcf file (default: *out*)
1. **-g|--annotation**  [FILE]: GTF file. Enables annotation of SV calls and fusion identification.
1. **-s|--min_support** [INT]: The minimum number of supporting reads required to be considered a SV (default: *2*)
1. **-m|--min_length** [INT]: Minimum SV length (default: *30*)
1. **-M|--max_length** [INT]: Maximum SV length (default: *20000*)

### Additional Parameters ###
1. **-h|--help**: Shows help message.
1. **-v|--version**: Shows current version.
1. **-p|--print_reads**: Print all contigs and associated reads as additional output *out.vcf.reads* ( *out* is the prefix specified in **-o|--output** )when activated.
1. **-H|--heuristic**: Use clustering heuristic when actived. Good for data with PCR duplicates.
1. **-D|--dump_contigs**: Dump contigs in fastq format for mapping.


You can always use 
```
./svict -h
```
to get a description of all the parameters. 

---

# Publication
See the [publication page](https://github.com/vpc-ccg/svict/blob/master/PUBLICATION.md) for citation and more information.


# Contact & Support
To report any bugs or issues please refer to the [issues page](https://github.com/vpc-ccg/svict/issues).
