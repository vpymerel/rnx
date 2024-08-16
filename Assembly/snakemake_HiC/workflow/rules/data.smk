
################################################################################
# sampling                                                                     #
################################################################################

rule sampling:
    input:
        'data/{sample}_{fr}.fastq.gz'
    output:
        'data/{sample}_{fr}_sampled.fastq.gz'
    log:
        'logs/{sample}_{fr}.sampling'
    conda:
        '../envs/sektq.yaml'
    threads:
        2
    resources:
        mem_mb = 8000
    params:
        runtime = '0-06:00:00'
    shell:
        "seqtk sample -s100 {input} {config[sampling_prop]} > data/{wildcards.sample}_{wildcards.fr}_sampled.fastq;"
        "gzip data/{wildcards.sample}_{wildcards.fr}_sampled.fastq;"

################################################################################
# links                                                                     #
################################################################################

rule links:
    input:
       'data/{sample}_{fr}.fastq.gz' if config["sampling_prop"]==1
       else 
       'data/{sample}_{fr}_sampled.fastq.gz'
    output:
        'results/{sample}/{sample}_{fr}_UnTrimmed.fastq.gz'
    log:
        'logs/{sample}_{fr}.links'
    threads:
        1
    resources:
        mem_mb = 1000
    params:
        runtime = '0-00:30:00'
    shell:
        "cp {input} {output};"
        
################################################################################
# links                                                                     #
################################################################################
        
rule trimming:
    input:
        expand('results/{{sample}}/{{sample}}_{fr}_UnTrimmed.fastq.gz',
        fr = FR)
    output:
        R1 = 'results/{sample}/{sample}_R1.fastq.gz',
        R2 = 'results/{sample}/{sample}_R2.fastq.gz'
    log:
        'logs/{sample}.trimming'
    conda:
        '../envs/trim_galore.yaml'
    threads:
        8
    resources:
        mem_mb = 16000
    params:
        runtime = '0-08:00:00'
    shell:
        "trim_galore -q 30 \
        --nextera \
        --length 100 \
        --cores 8 \
        --fastqc \
        --output_dir results/{wildcards.sample}/trimming \
        --paired {input};"
        "mv results/{wildcards.sample}/trimming/{wildcards.sample}_R1_UnTrimmed_val_1.fq.gz {output.R1};"
        "mv results/{wildcards.sample}/trimming/{wildcards.sample}_R2_UnTrimmed_val_2.fq.gz {output.R2};"

        