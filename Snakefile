configfile: "config.yml"

rule gen_reads:
    params:
        name=config['NAME'],
        ref=config['REF'],
        lengths=config['LENGTHS'],
        n_seqs=config['N_SEQS'],
        n_workers=config['N_WORKERS'],
        output_dir="simulations"
    threads:
        95
    shell:
        "mkdir -p {params.output_dir} && "
        "python sim.py --name {params.name} --ref {params.ref} --lengths {params.lengths} --n_seqs {params.n_seqs} --n_workers {params.n_workers} --out {params.output_dir}"

# rule all:
    # input:
        

rule bwa_map:
    input:
        fa="/storage/data/resources/reference_genome/hg38/BWAIndex/genome.fa",
        fastqs=lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 1
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.fa} {input.fastqs} |"
        "samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.srt.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule merge_bams:
    input:
        "sorted_reads/{sample}.srt.bam"
    output:
        bam="merged.bam"
    threads: 90
    shell:
        "samtools merge -@ {threads} {output.bam} {wildcards.sample}"
        "samtools index -@ {threads} {output.bam}"
