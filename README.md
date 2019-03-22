**SViCT**: Structural Variant detrction in Circulating Tumor DNA
===================
### What is SViCT?
SViCT is a computational tool for detecting structural variations from targeted sequencing reads derived from cell free DNA (cfDNA) containing low dilutions of circulating tumor DNA (ctDNA).

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
