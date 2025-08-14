from snakemake.io import glob_wildcards
SAMPLES = glob_wildcards("results/trimmed/{sample}_1P.fastq.gz").sample
REFERENCE = "reference/NC_045512.2.fasta"

SAMPLES = ["SRR23609077","SRR23609078","SRR23609079","SRR23609080","SRR23609081","SRR23609082","SRR23609083","SRR23609084","SRR23609085","SRR23609086"]
REF = "reference/NC_045512.2.fasta"

rule all:
    input:
        expand("results/fastqc_trimmed/{sample}_1P_fastqc.zip", sample=SAMPLES),
        expand("results/fastqc_trimmed/{sample}_2P_fastqc.zip", sample=SAMPLES),
        "results/multiqc/multiqc_report.html",
        expand("results/bam/{sample}.dedup.bam", sample=SAMPLES),
        expand("results/bam/{sample}.dedup.bam.bai", sample=SAMPLES),
        expand("results/bam/{sample}.markdup.txt", sample=SAMPLES),
        expand("results/coverage/{sample}.depth.txt", sample=SAMPLES),
        expand("results/vcf/{sample}.norm.vcf.gz", sample=SAMPLES),
        expand("results/vcf/{sample}.norm.vcf.gz.tbi", sample=SAMPLES),
        "results/vcf/merged.vcf.gz",
        "results/vcf/merged.vcf.gz.tbi"

rule faidx:
    input:
        REF
    output:
        REF + ".fai"
    shell:
        "samtools faidx {input}"

rule bwa_index:
    input:
        REF
    output:
        REF + ".amb",
        REF + ".ann",
        REF + ".bwt",
        REF + ".pac",
        REF + ".sa"
    shell:
        "bwa index {input}"

rule fastqc_trimmed:
    input:
        r1="results/trimmed/{sample}_1P.fastq.gz",
        r2="results/trimmed/{sample}_2P.fastq.gz"
    output:
        "results/fastqc_trimmed/{sample}_1P_fastqc.html",
        "results/fastqc_trimmed/{sample}_1P_fastqc.zip",
        "results/fastqc_trimmed/{sample}_2P_fastqc.html",
        "results/fastqc_trimmed/{sample}_2P_fastqc.zip"
    shell:
        "mkdir -p results/fastqc_trimmed && fastqc -q -o results/fastqc_trimmed {input.r1} {input.r2}"

rule multiqc:
    input:
        expand("results/fastqc_trimmed/{sample}_1P_fastqc.zip", sample=SAMPLES),
        expand("results/fastqc_trimmed/{sample}_2P_fastqc.zip", sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        "mkdir -p results/multiqc && multiqc -q -o results/multiqc results/fastqc_trimmed"

rule map_sort:
    input:
        ref=REF,
        fai=REF + ".fai",
        amb=REF + ".amb",
        ann=REF + ".ann",
        bwt=REF + ".bwt",
        pac=REF + ".pac",
        sa=REF + ".sa",
        r1="results/trimmed/{sample}_1P.fastq.gz",
        r2="results/trimmed/{sample}_2P.fastq.gz"
    output:
        "results/bam/{sample}.sorted.bam"
    threads: 6
    shell:
        "mkdir -p results/bam && "
        "bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' {input.ref} {input.r1} {input.r2} | "
        "samtools sort -@ {threads} -o {output} -"

rule index_bam:
    input:
        "results/bam/{sample}.sorted.bam"
    output:
        "results/bam/{sample}.sorted.bam.bai"
    threads: 2
    shell:
        "samtools index -@ {threads} {input}"

rule name_sort:
    input:
        "results/bam/{sample}.sorted.bam"
    output:
        "results/bam/{sample}.namesort.bam"
    threads: 4
    shell:
        "samtools sort -n -@ {threads} -o {output} {input}"

rule fixmate:
    input:
        "results/bam/{sample}.namesort.bam"
    output:
        "results/bam/{sample}.fixmate.bam"
    threads: 2
    shell:
        "samtools fixmate -@ {threads} -m {input} {output}"

rule pos_sort_fixmate:
    input:
        "results/bam/{sample}.fixmate.bam"
    output:
        "results/bam/{sample}.fixmate.pos.bam"
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule markdup:
    input:
        "results/bam/{sample}.fixmate.pos.bam"
    output:
        bam="results/bam/{sample}.dedup.bam",
        bai="results/bam/{sample}.dedup.bam.bai",
        metrics="results/bam/{sample}.markdup.txt"
    threads: 4
    params:
        tmp="results/bam/tmp/{sample}"
    shell:
        "mkdir -p {params.tmp} && "
        "samtools markdup -@ {threads} -T {params.tmp}/md {input} {output.bam} 2> {output.metrics} && "
        "samtools index -@ {threads} {output.bam}"

rule coverage_depth:
    input:
        "results/bam/{sample}.dedup.bam"
    output:
        "results/coverage/{sample}.depth.txt"
    shell:
        "mkdir -p results/coverage && samtools depth -aa -d 0 {input} > {output}"

rule bcftools_call:
    input:
        bam="results/bam/{sample}.dedup.bam",
        bai="results/bam/{sample}.dedup.bam.bai",
        ref=REF,
        fai=REF + ".fai"
    output:
        vcf="results/vcf/{sample}.vcf.gz",
        tbi="results/vcf/{sample}.vcf.gz.tbi"
    threads: 4
    shell:
        "mkdir -p results/vcf && "
        "bcftools mpileup -f {input.ref} {input.bam} -Ou | "
        "bcftools call -mv -Oz -o {output.vcf} && "
        "bcftools index -t {output.vcf}"

rule bcftools_norm:
    input:
        vcf="results/vcf/{sample}.vcf.gz",
        ref=REF
    output:
        vcf="results/vcf/{sample}.norm.vcf.gz",
        tbi="results/vcf/{sample}.norm.vcf.gz.tbi"
    threads: 2
    shell:
        "bcftools norm -f {input.ref} -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index -t {output.vcf}"

rule merge_vcf:
    input:
        expand("results/vcf/{sample}.norm.vcf.gz", sample=SAMPLES)
    output:
        vcf="results/vcf/merged.vcf.gz",
        tbi="results/vcf/merged.vcf.gz.tbi"
    shell:
        "bcftools merge -Oz -o {output.vcf} {input} && bcftools index -t {output.vcf}"

