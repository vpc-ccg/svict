**SVICT**: Structural Variant detrction in Circulating Tumor DNA
===================
### What is SVICT?
SVICT is a computational tool for detecting structural variations from targeted sequencing reads derived from cell free DNA (cfDNA) containing low dilutions of circulating tumor DNA (ctDNA).

### How do I get SVICT?
Just clone our repository and issue the `make` command:
```
git clone --recursive https://github.com/vpc-ccg/cfdna-sv.git
cd cfdna-sv && make
```

> **Note**: You will need at least g++ 4.9 to compile the source code.

### How do I run SVICT?
You can use 
```
cfdna-sv/SVICT -h
```
to get a description of each parameter. 

#### Simulation Datasets Used for Evaluation
Please check [this link](https://goo.gl/PTzJec) to download the simulation datasets that we used for evaluating SVICT. The folder contains 3 data files:
1. sim.pe.sorted.bam: BAM file containing simulated cfDNA reads from a Venter genome with inserted SVs
2. sim.pe.sorted.bam.bai: Corresponding index
3. Homo_sapiens.GRCh38.87.dna.chromosomes.fa: reference genome for chromosome 1 from GRCh38.

We also provide a checksum file ***md5.sum*** for checking file integrity.

To run SVICT on this dataset, type

```
./SVICT -i sim.pe.sorted.bam -r Homo_sapiens.GRCh38.87.dna.chromosomes.fa -o [output_prefix]
```

Feel free to change **output_prefix** to any name you want. Remember to specify the paths correctly for ***Homo_sapiens.GRCh38.87.dna.chromosomes.fa*** and ***sim.pe.sorted.bam*** if you download these files to a folder other than the one you run SVICT. The VCF file with the prediction results will be generated as **output_prefix.vcf** in the current directory.


---


### Contact & Support

Feel free to drop any inquiry to [agawrons at sfu dot ca](mailto:).
