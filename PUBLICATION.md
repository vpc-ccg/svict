**SViCT**: Structural Variant detrction in Circulating Tumor DNA
===================
SViCT is a computational tool for detecting structural variations from cell free DNA (cfDNA) containing low dilutions of circulating tumor DNA (ctDNA).

# Installation 

## BIOCONDA

## Installation from Source


To install SViCT, you first need to fetch it from the github repository. After downloading, change the current directory to the source directory 'svict' and run make in the terminal to create the binary file 'svict'.

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

You can also go to releases page, click on the desired version and then click either **Source Code (zip)** or **Source code (tar.gz)** link to download the file. After decompressing it, you just switch to the ```svict``` directory and run ```make```.

### Prerequisite ###
You will need at least g++ 4.9 to compile the source code.


# Running SViCT
SViCT requires **coordinate-sorted BAM/SAM** and **reference genome FASTA** files to run detections:

```
./svict -i [Sorted BAM/SAM] -r [Reference Genome FASTA]
```
By default, the output is written to "out.vcf" in the current folder.

## Test Datasets
To grab sample data and test ```SViCT```, please download the following two files:
```
curl -L https://ndownloader.figshare.com/files/12380225 --output sim.75.sorted.bam
curl -L https://ndownloader.figshare.com/files/10144653 --output Homo_sapiens.GRCh38.87.dna.chromosomes.fa
```

Run the following commands
```
./svict -i sim.75.sorted.bam -r Homo_sapiens.GRCh38.87.dna.chromosomes.fa -o sim.75.vcf
```
The VCF file with the prediction results will be generated as **sim.75.vcf** in the current directory, and it shoule be identical to ```out.75.vcf```. 

> You need to the paths correctly for ***Homo_sapiens.GRCh38.87.dna.chromosomes.fa*** and ***sim.100.sorted.bam*** if you download these files to a folder other the ```svict``` folder. 


## Command Options ## 
### Mandatory Parameters ###
1. **-i|--input**: Input BAM/SAM file sorted by coordinate
1. **-r|--reference**: Reference genome fasta file which was used for the above mapping file

### Main Optional Parameters ###
1. **-o|--output**: Prefix of output vcf file (default: *out*)
1. **-g|--annotation**: GTF file which enables annotation of SV calls and fusion identification
1. **-s|--min_support**: The minimum number of supporting reads required to be considered a SV (default: *2*)
1. **-S|--max_support**: The maximum number of supporting reads required to be considered a SV (default: unlimited)
1. **-m|--min_length**: Minimum SV length (default: *30*)
1. **-M|--max_length**: Maximum SV length (default: *20000*)

### Additional Parameters ###
1. **-p|--print_reads**: Print all contigs and associated reads as additional output.
1. **-P|--print_stats**: Print statistics to stderr.
1. **-w|--window_size**: Clustering window (default: *3*).
1. **-d|--min_sc**: Minimum soft clip to consider (default: *10*).
1. **-n|--no_indel**: Disable indel parsing (```I``` and ```D``` in cigar).
1. **-O|--assembler_overlap**: Required read overlap for assembly (default: *50*).
1. **-a|--anchor**: Anchor length (default: *30*).
1. **-k|--kmer**: k-mer length (default: *14*).
1. **-u|--uncertainty**: Uncertainty (default: *8*).
1. **-c|--sub_optimal**: Maximum difference from longest path (default: *0* - co-optimals only, negative value disables).
1. **-H|--heuristic**: Use clustering heuristic when actived. Good for data with PCR duplicates.
1. **-D|--dump_contigs**: Dump contigs in fastq format for mapping.
1. **-R|--resume**: Resume at the interval chaining stage with mapped contigs.


### How do I get SViCT?
Just clone our repository and issue the `make` command:
```
git clone https://github.com/vpc-ccg/svict.git
cd svict && make
```

> **Note**: You will need at least g++ 4.9 to compile the source code.

### How do I run SViCT?

At minimum, SViCT requires two parameters
```
./svict -i [Sorted BAM/SAM] -r [Reference Genome FASTA]
```
and output is written to "out.vcf".

You can use 
```
svict/svict -h
```
to get a description of all the parameters. 

#### Simulation Datasets Used for Evaluation
Please check [this link]( https://figshare.com/articles/Simulation_Datasets_for_Evaluation/5758539 ) to download the simulation datasets that we used for evaluating SViCT. The folder contains 4 data files:
1. sim.150.sorted.bam: BAM file containing simulated 2*150bp cfDNA reads from a Venter genome with inserted SVs
2. sim.100.sorted.bam: BAM file containing simulated 2*100bp cfDNA reads from a Venter genome with inserted SVs
3. sim.75.sorted.bam: BAM file containing simulated 2*75bp cfDNA reads from a Venter genome with inserted SVs
4. Homo_sapiens.GRCh38.87.dna.chromosomes.fa: GRCh38 reference genome for the above mappings.

We also provide a checksum file ***md5.sum*** for checking file integrity.

To run SViCT on this dataset, type

```
./svict -i sim.100.sorted.bam -r Homo_sapiens.GRCh38.87.dna.chromosomes.fa
```

Remember to specify the paths correctly for ***Homo_sapiens.GRCh38.87.dna.chromosomes.fa*** and ***sim.100.sorted.bam*** if you download these files to a folder other than the one you run SViCT. The VCF file with the prediction results will be generated as **out.vcf** in the current directory. An alternate output prefix can be specified with "-o".


### Existing SV callers Used in Performance Comparisons ###
1. [Lumpy2](https://github.com/arq5x/lumpy-sv)   0.2.13 
2. [GRIDSS](https://github.com/PapenfussLab/gridss) 1.4.3
3. [Socrates](https://github.com/jibsch/Socrates) 1.13.1
4. [Pindel](https://github.com/genome/pindel) 0.2.5
5. [Delly2](https://github.com/dellytools/delly) 0.7.8
---

### Publications
**Structural variation and fusion detection using targeted sequencing data from circulating cell free DNA.** Alexander R Gawroński, Yen-Yi Lin, Brian McConeghy,   Stephane LeBihan, Hossein Asghari, Can Koçkan, Baraa Orabi, Nabil Adra, Roberto Pili, Colin C Collins, S Cenk Sahinalp, Faraz Hach. [Nucleic Acids Research, 2019](https://doi.org/10.1093/nar/gkz067)


### Contact & Support

Feel free to drop any inquiry to [agawrons at sfu dot ca](mailto:).
