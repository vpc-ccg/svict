**SViCT**: Structural Variant detrction in Circulating Tumor DNA
===================
### What is SViCT?
SViCT is a computational tool for detecting structural variations from targeted sequencing reads derived from cell free DNA (cfDNA) containing low dilutions of circulating tumor DNA (ctDNA).

### How do I get SViCT?
Just clone our repository and issue the `make` command:
```
git clone https://github.com/vpc-ccg/cfdna-sv.git
cd cfdna-sv && make
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
cfdna-sv/svict -h
```
to get a description of all the parameters. 

#### Simulation Datasets Used for Evaluation
Please check [this link](https://goo.gl/PTzJec) to download the simulation datasets that we used for evaluating SViCT. The folder contains 2 data files:
1. sim.100.sorted.bam: BAM file containing simulated cfDNA reads from a Venter genome with inserted SVs (75bp and 150bp read data is also availible)
2. Homo_sapiens.GRCh38.87.dna.chromosomes.fa: reference genome from GRCh38.

We also provide a checksum file ***md5.sum*** for checking file integrity.

To run SViCT on this dataset, type

```
./svict -i sim.100.sorted.bam -r Homo_sapiens.GRCh38.87.dna.chromosomes.fa
```

Remember to specify the paths correctly for ***Homo_sapiens.GRCh38.87.dna.chromosomes.fa*** and ***sim.100.sorted.bam*** if you download these files to a folder other than the one you run SViCT. The VCF file with the prediction results will be generated as **out.vcf** in the current directory. An alternate output prefix can be specified with "-o".


---


### Contact & Support

Feel free to drop any inquiry to [agawrons at sfu dot ca](mailto:).
