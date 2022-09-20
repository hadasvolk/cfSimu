# cfSimu-wgsim
Cell free simulation using wgsim and the python modified version [pywgsim](https://github.com/ialbert/pywgsim)

## About
The python script, sim.py, generates paired-end fastq files based on hoobari's fetal and maternal fragment length distrubtion. The output is multiple fastq files depending on user defined number of sequences and number of cores. Please refer to script for input options.
Alignment and merging of bams can be achevied using the attached Snakefile

## Run
1. copy hoobari's output maternal_lengths.csv and fetal_lengths.csv to cwd
2. generate compressed fastqs with sim.py
```python 
python /shared/home/hadas/identifai/repos/cfSimu-wgsim/sim.py --name fam50_maternal --n_seqs 1000000000 -nw 95 -o reads_maternal
```
3. change Snakefile accordingly </br>
run </br>
```python
snakemae -c <n cores>
```