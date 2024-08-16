################################################################################
# fastqc                                                                       #
################################################################################


rule fastqc:
    input:
        R1 = 'results/{sample}/{sample}_R1.fastq.gz',
        R2 = 'results/{sample}/{sample}_R2.fastq.gz'
    output:
        directory('results/{sample}/fastqc/')
    log:
        'logs/{sample}.fastqc'
    conda:
        '../envs/fastqc.yaml'
    threads:
        8
    resources:
        mem_mb = 8000
    params:
        runtime = '0-01:00:00'
    shell:
        "mkdir -p {output};"
        "fastqc {input.R1} {input.R2} -o {output} -t {threads} --extract;"


################################################################################
# bamqc                                                             #
################################################################################

rule bam_sorting:
    input:
        'results/{sample}/mapped/{sample}_{fr}.bam'
    output:
        'results/{sample}/mapped/{sample}_{fr}.sorted.bam'
    log:
        'logs/{sample}_{fr}.bam_sorting'
    conda:
        '../envs/samtools.yaml'
    threads:
        8
    resources:
        mem_mb = 16000
    params:
        runtime = '0-06:00:00'
    shell:
        "samtools sort -@ 8 {input} -o {output};"

rule bamqc:
    input:
        'results/{sample}/mapped/{sample}_{fr}.sorted.bam'
    output:
        'results/{sample}/mapped/{sample}_{fr}.sorted_stats/qualimap.pdf'
    log:
        'logs/{sample}_{fr}.bamqc'
    threads:
        8
    resources:
        mem_mb = 16000
    params:
        runtime = '0-01:00:00'
    shell:
        "rm -rf results/{wildcards.sample}/mapped/{wildcards.sample}_{wildcards.fr}.sorted_stats;"
        "module load gcc/10.4.0 qualimap/2.2.1;"
        "qualimap bamqc -nt 8 -bam {input} -outfile qualimap.pdf"
