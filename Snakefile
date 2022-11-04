import os
import pandas as pd

SAMPLE=config["SAMPLE"]
NAME=config["NAME"]
OUTPUT=config["OUTPUT"]
N_SEQS=config["N_SEQS"]

OUTDIR="{}_sims".format(SAMPLE)
REF="genome.fa"
SIM="sim.py"
CONDA="cfsimu.yml"
TEMP="/tmp"
THREADS=45

# Create output directory
try:
    os.path.isfile("{}_lengths.csv".format(SAMPLE))
except IOError:
    print("{}_lengths.csv missing in cwd".format(SAMPLE))
    os.exit(1)

lengths_file = "{}_lengths.csv".format(SAMPLE)

rule all:
    input:
        os.path.join(OUTPUT, SAMPLE, "{}.srt.cram".format(SAMPLE)),
        os.path.join(OUTPUT, SAMPLE,"{}.srt.cram.stats".format(SAMPLE)),
        os.path.join(OUTPUT, SAMPLE,"{}.mosdepth.summary.txt".format(SAMPLE)),
        # os.path.join(OUTPUT, SAMPLE,"done.txt")

rule simulator:
    params:
        name=NAME,
        n_seqs=N_SEQS,
        sim=SIM,
        out=OUTDIR,
        lengths=lengths_file
    # conda:
    #     CONDA
    output:
        os.path.join(OUTPUT, OUTDIR, "done.txt")
    threads:
        THREADS
    resources:
        tmpdir=TEMP
    shell:
        "python {params.sim} --name {params.name} --n_seqs {params.n_seqs} -nw {threads} -o {params.out} -l {params.lengths}"

rule merge_fastqs_r1:
    input:
        os.path.join(OUTPUT, OUTDIR, "done.txt")
    output:
        temp(os.path.join(OUTPUT, SAMPLE, "{}_R1.fastq.gz".format(SAMPLE)))
    params:
        out=OUTDIR
    shell:
        "cat {params.out}/read1.*.gz > {output}"

rule merge_fastqs_r2:
    input:
        os.path.join(OUTPUT, OUTDIR, "done.txt")
    output:
        temp(os.path.join(OUTPUT, SAMPLE, "{}_R2.fastq.gz".format(SAMPLE)))
    params:
        out=OUTDIR
    shell:
        "cat {params.out}/read2.*.gz > {output}"

rule bwa_map:
    input:
       fa=REF,
       r1=os.path.join(OUTPUT, SAMPLE, "{}_R1.fastq.gz".format(SAMPLE)),
       r2=os.path.join(OUTPUT, SAMPLE, "{}_R2.fastq.gz".format(SAMPLE))
    output:
        os.path.join(OUTPUT, SAMPLE,"{}.srt.cram".format(SAMPLE))
    params:
        rg=r"@RG\tID:{}\tSM:{}".format(SAMPLE,SAMPLE)
    log:
        os.path.join(OUTPUT, SAMPLE,"logs/bwa_mem/{}.log".format(SAMPLE))
    threads: 
        THREADS
    resources:
        tmpdir=TEMP
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.fa} {input.r1} {input.r2} |"
        "samtools sort -O bam -l 0 - |"
        "samtools view -T {input.fa} -C -o {output} -) 2> {log};"
        "samtools index -@ {threads} {output}"

rule samtools_stats:
    input:
        os.path.join(OUTPUT, SAMPLE,"{}.srt.cram".format(SAMPLE))
    output:
        os.path.join(OUTPUT, SAMPLE,"{}.srt.cram.stats".format(SAMPLE))
    threads: 
        THREADS
    resources:
        ref=REF
    shell:
        "samtools stats -@ {threads} --reference {resources.ref} {input} > {output}"

rule mosdepth:
    input:
        os.path.join(OUTPUT, SAMPLE,"{}.srt.cram".format(SAMPLE))
    output:
        os.path.join(OUTPUT, SAMPLE,"{}.mosdepth.summary.txt".format(SAMPLE))
    params:
        ref=REF,
        prefix=os.path.join(OUTPUT, SAMPLE,"{}".format(SAMPLE))
    threads: 
        THREADS
    shell:
        "mosdepth -n -f {params.ref} -t {threads} {params.prefix} {input}"

rule cleanup:
    input:
        os.path.join(OUTPUT, SAMPLE,"{}.srt.cram".format(SAMPLE))
    output:
        os.path.join(OUTPUT, SAMPLE, "done.txt")
    params:
        out=OUTDIR
    shell:
        "rm -rf {params.out};"
        "touch {output}"
