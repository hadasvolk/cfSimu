import os
import pandas as pd

SAMPLE='fetal'
OUTPUT='/storage/users/hadas/chromo/wilchrom_sims/sim04'
FASTQS=os.path.join(OUTPUT, "fetal_reads/fastq_dict.tsv")

REF="/storage/data/resources/reference_genome/hg38/BWAIndex/genome.fa"

samples_table = pd.read_csv(FASTQS, sep="\t").set_index("name", drop=False)

# fastq1 input function definition
def fq1_from_sample(wildcards):
    return samples_table.loc[wildcards.sample, "read1"]

# fastq2 input function definition
def fq2_from_sample(wildcards):
    return samples_table.loc[wildcards.sample, "read2"]

rule all:
    input:
        expand(os.path.join(OUTPUT, SAMPLE,"mapped_reads/{sample}.cram"), sample=samples_table.name),
        os.path.join(OUTPUT, SAMPLE, "fetal.srt.cram")

rule bwa_map:
    input:
       fa=REF,
       r1=fq1_from_sample,
       r2=fq2_from_sample
    output:
        temp(os.path.join(OUTPUT, SAMPLE,"mapped_reads/{sample}.cram"))
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        os.path.join(OUTPUT, SAMPLE,"logs/bwa_mem/{sample}.log")
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.fa} {input.r1} {input.r2} |"
        "samtools sort -O bam -l 0 -T /tmp - |"
        "samtools view -T {input.fa} -C -o {output} -) 2> {log}"

rule merge_bams:
    output:
        os.path.join(OUTPUT, SAMPLE, "fetal.srt.cram")
    params:
        d=os.path.join(OUTPUT, SAMPLE,"mapped_reads/*.cram")
    threads: 90
    shell:
        "samtools merge -@ {threads} {output} {params.d};"
        "samtools index -@ {threads} {output};"
        "samtools stats {output} > {output}.stats"
