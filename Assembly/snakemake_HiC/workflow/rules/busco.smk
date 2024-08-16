################################################################################
# busco                                                                        #
################################################################################

rule busco:
    input:
        'results/{sample}/{method}.fa'
    output:
        'results/{sample}/{method}/busco/run_insecta_odb10/full_table.tsv'
    log:
        'logs/{sample}.{method}.busco'
    container:
        "docker://ezlabgva/busco:v5.4.5_cv1"
    threads:
        12
    resources:
        mem_mb = 24000
    params:
        runtime = '0-12:00:00'
    shell:
        "cd results/{wildcards.sample}/{wildcards.method}/;"
        "busco -m genome \
        -i ../../{input} \
        -o busco \
        -l insecta_odb10 \
        -f \
        -c {threads};"
        
################################################################################
# backmapping                                                                  #
################################################################################
                
rule backmapping:
    input:
        reads = 'data/{sample}.fastq.gz',
        assembly = 'results/{sample}/{method}.fa'
    output:
        bam = 'results/{sample}/{method}/backmapping/mapped.bam',
        bai = 'results/{sample}/{method}/backmapping/mapped.bam.bai'
    log:
        'logs/{sample}.{method}.backmapping'
    conda:
        '../envs/purge_haplotigs.yaml'
    threads:
        12
    resources:
        mem_mb = 96000
    params:
        ali_tmp = temp('results/{sample}/{assembler}/backmapping/{step}_tmp.ali'),
        runtime = '0-06:00:00'
    shell:
        "minimap2 -t {threads} \
        -ax map-hifi \
        {input.assembly} \
        {input.reads} \
        --secondary=no | samtools sort -m 24G \
        -o {output.bam} \
        -T {params.ali_tmp};"
        "samtools index {output.bam};"