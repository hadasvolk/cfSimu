# cfSimu
CfSimu is the first of its kind cell free fragment length distribution estimator and simulator. The tool Allows for truth like profiling of the fragment lengths distribution of paired-end NGS sequencing fragments using least mean square fitting for range of non-linear functions. The tool performs fit assessment and generate simulated paired end short read sequences in fastq and aligned formats. 

Cell free simulation engine wgsim and the python modified version [pywgsim](https://github.com/ialbert/pywgsim)

# UNDER HEAVY DEVELOPMENT
## About
The python script, sim.py, generates paired-end fastq files based a user given fetal and maternal fragment length distrubtion as a tabular csv file. The output is multiple fastq files depending on user defined number of sequences and number of cores. Please refer to script for input options (python sim.py -h)

Alignment and merging of crams can be achevied using the attached Snakefile

Generate maternal and fetal alignmets based on fetal fraction prefrences and merge with samtools 

## Run
1. Create a csv read distrubtion file with the following columns: 'length', 'COUNT(length)'
2. generate compressed fastqs with sim.py
```python 
python src/sim.py --name maternal_family --n_seqs 1000000000 -nw 95 -o reads_maternal
```
3. change Snakefile accordingly </br>
run </br>
```python
snakemae -c <n cores>
```
