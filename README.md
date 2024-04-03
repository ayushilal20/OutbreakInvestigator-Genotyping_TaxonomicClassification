# Pipeline for genotyping and taxonomy: 

## Tools used for the pipeline: 

1. Geneotyping: 
For genotyping we used MAAST and MLST

"This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust".​

Seemann T, mlst Github https://github.com/tseemann/mlst 

Shi, Z.J., Nayfach, S. & Pollard, K.S. Maast: genotyping thousands of microbial strains efficiently. Genome Biol 24, 186 (2023). https://doi.org/10.1186/s13059-023-03030-8

2. Taxonomic Analysis: 
For taxonomy we used kraken2, lissero, and blastn

Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0

Josh Zhang, LisSero Github https://github.com/MDU-PHL/LisSero

Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+: architecture and applications. BMC Bioinformatics 10, 421 (2009). https://doi.org/10.1186/1471-2105-10-421

3. Genome Quality Assessment

Parks, Donovan H., et al. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome research 25.7 (2015): 1043-1055.

Gurevich A, Saveliev V, Vyahhi N, Tesler G. QUAST: quality assessment tool for genome assemblies. Bioinformatics. 2013 May 1;29(8):1072-5. doi: 10.1093/bioinformatics/btt086. Epub 2013 Feb 21. PMID: 23426934.


## Pipeline:
The pipeline takes the input as the fasta files from the assembly and those are used by all programs. Kraken2 and LisSero outputs are used to determine which reference genome to use for Quast 

## Pipeline Usage: 

```python
sh project3_script.sh [input_dir]
```

The pipeline fully automates the cumbersome process and it generates the required files that will help in taxonomic classification and comparative genomics.  
