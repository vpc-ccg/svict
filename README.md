**SViCT**: Structural Variant detection in Circulating Tumor DNA
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
conda install -c bioconda svict
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
SViCT requires **coordinate-sorted BAM/SAM** generated from Illumina reads (short read technologies) and **reference genome FASTA** files to run detections:

```
./svict -i [input] -r [reference]
```

The output will be written to **out.vcf** in the current folder.

## Running Example
To test ```SViCT```, please first download the following two files ( ~ 3GB in total):
```
curl -L https://ndownloader.figshare.com/files/12380225 --output sim.75.sorted.bam
curl -L https://ndownloader.figshare.com/files/10144653 --output Homo_sapiens.GRCh38.87.dna.chromosomes.fa
```

Type the following command to run ```svict```:
```
./svict -i sim.75.sorted.bam -r Homo_sapiens.GRCh38.87.dna.chromosomes.fa -o out
```
The VCF file with the prediction results will be generated as **out.vcf** in the current directory. The output should be the same as [**original result**](https://ndownloader.figshare.com/files/14677538).


You can always use 
```
./svict -h
```
to get a description of all the parameters. 

---

# Publication
If you use SViCT in your publication, please cite the following article:

**Structural variation and fusion detection using targeted sequencing data from circulating cell free DNA.** Alexander R Gawroński, Yen-Yi Lin, Brian McConeghy,   Stephane LeBihan, Hossein Asghari, Can Koçkan, Baraa Orabi, Nabil Adra, Roberto Pili, Colin C Collins, S Cenk Sahinalp, Faraz Hach. [Nucleic Acids Res. 2019 Apr 23;47(7):e38. doi: 10.1093/nar/gkz067](https://doi.org/10.1093/nar/gkz067)


See the [publication page](https://github.com/vpc-ccg/svict/blob/master/PUBLICATION.md) for details about the experiements.


# Contact & Support
To report any bugs or issues please refer to the [issues page](https://github.com/vpc-ccg/svict/issues).
