# cfSimu-wgsim
Cell free simulation using Wgsim

# Run
1. copy hoobari's output maternal_lengths.csv and fetal_lengths.csv to cwd
2. generate compressed fastqs with sim.py
'''python 
python /shared/home/hadas/identifai/repos/cfSimu-wgsim/sim.py --name fam50_maternal --n_seqs 1000000000 -nw 95 -o reads_maternal
'''
3. change Snakefile accordingly
run
'''python
snakemae -c <n cores>
'''