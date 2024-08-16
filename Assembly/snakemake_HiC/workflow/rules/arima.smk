# Arima style mapping                                                          #


################################################################################
# bwa index                                                        #
################################################################################

rule bwa_indexing:
    input:
        'data/{sample}.fasta'
    output:
        'results/{sample}/mapped/{sample}.fasta',
    log:
        'logs/{sample}.bwa_indexing'
    conda:
        '../envs/bwa.yaml'
    threads:
        4
    resources:
        mem_mb = 80000
    params:
        runtime = '0-01:00:00',
        dir = 'results/{sample}/mapped/'
    shell:
        "rm -rf {params.dir}; mkdir -p {params.dir};"
        "cp {input} {params.dir}/{wildcards.sample}.fasta;"
        "cd {params.dir};"
        "bwa index {wildcards.sample}.fasta;"

################################################################################
# samtools index                                                        #
################################################################################

rule samtools_indexing:
    input:
        'results/{sample}/mapped/{sample}.fasta',
    output:
        'results/{sample}/mapped/{sample}.fasta.fai'
    log:
        'logs/{sample}.samtools_indexing'
    conda:
        '../envs/samtools.yaml'
    threads:
        4
    resources:
        mem_mb = 80000
    params:
        runtime = '0-01:00:00'
    shell:
        "samtools faidx {input};"
        
################################################################################
# bwa mem                                                   #
################################################################################

rule mapping:
    input:
        assembly = 'results/{sample}/mapped/{sample}.fasta',
        reads = 'results/{sample}/{sample}_{fr}.fastq.gz'
    output:
        bam = 'results/{sample}/mapped/{sample}_{fr}.bam'
    log:
        'logs/{sample}_{fr}.mapping'
    conda:
        '../envs/mapping.yaml'
    threads:
        16
    resources:
        mem_mb = 32000
    params:
        runtime = '0-24:00:00',
    shell:
        "bwa mem -t {threads} {input.assembly} {input.reads} | samtools view -@ {threads} -Sb - > {output};"

################################################################################
# filtering chimerics                                                   #
################################################################################
    
rule filter_chimerics:
    input:
        'results/{sample}/mapped/{sample}_{fr}.bam'
    output:
        'results/{sample}/mapped/{sample}_{fr}_filtered.bam'
    log:
        'logs/{sample}_{fr}.filter_chimerics'
    conda:
        '../envs/filter_chimerics.yaml'
    threads:
        4
    resources:
        mem_mb = 8000
    params:
        runtime = '0-12:00:00',
    shell:
        "samtools view -h {input} | perl resources/filter_five_end.pl | samtools view -Sb - > {output};"

################################################################################
# combining                                                    #
################################################################################


rule combining:
    input:
        faidx = 'results/{sample}/mapped/{sample}.fasta.fai',
        bams = expand('results/{sample}/mapped/{sample}_{fr}_filtered.bam',
        sample = SAMPLES,
        fr = FR)
    output:
        'results/{sample}/mapped/{sample}.bam'
    log:
        'logs/{sample}.combining'
    conda:
        '../envs/combining.yaml'
    threads:
        16
    resources:
        mem_mb = 100000
    params:
        runtime = '0-24:00:00',
        dir = 'results/{sample}/mapped/'
    shell:
        "perl resources/two_read_bam_combiner.pl \
        {input.bams} \
        samtools \
        {config[QC]} | samtools view \
        -bS \
        -t {input.faidx} - | samtools sort \
        -@ {threads} \
        -o {output} -"

################################################################################
# add_RG                                                    #
################################################################################

rule add_RG:
    input:
        'results/{sample}/mapped/{sample}.bam'
    output:
        'results/{sample}/mapped/{sample}_rg.bam'
    log:
        'logs/{sample}.add_RG'
    conda:
        '../envs/picard.yaml'
    threads:
        4
    resources:
        mem_mb = 8000
    params:
        runtime = '0-24:00:00',
    shell:
        "mkdir -p temp/{wildcards.sample};"
        "java -Xmx4G -Djava.io.tmpdir=temp/{wildcards.sample} -jar resources/picard.jar AddOrReplaceReadGroups \
        INPUT={input} \
        OUTPUT={output} \
        ID={wildcards.sample} \
        LB={wildcards.sample} \
        SM={wildcards.sample} \
        PL=ILLUMINA PU=none"

################################################################################
# mark_duplicates                                                    #
################################################################################

rule mark_duplicates:
    input:
        'results/{sample}/mapped/{sample}_rg.bam'
    output:
        bam = 'results/{sample}/{sample}.bam',
        metrics = 'results/{sample}/mapped/{sample}.metrics.txt'
    log:
        'logs/{sample}.mark_duplicates'
    conda:
        '../envs/picard.yaml'
    threads:
        8
    resources:
        mem_mb = 40000
    params:
        runtime = '0-24:00:00',
        dir = 'results/{sample}/mapped/'
    shell:
        "mkdir -p temp/{wildcards.sample};"
        "mkdir -p plop/{wildcards.sample};"
        "java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/{wildcards.sample} -jar resources/picard.jar MarkDuplicates \
        INPUT={input} \
        OUTPUT={output.bam} \
        METRICS_FILE={output.metrics} \
        TMP_DIR=plop/{wildcards.sample} \
        ASSUME_SORTED=TRUE \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=TRUE;"