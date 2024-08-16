################################################################################
# salsa2_prep                                                    #
################################################################################

rule salsa2_prep:
    input:
        assembly = 'data/{sample}.fasta',
        aln = 'results/{sample}/{sample}.bam'
    output:
        dir = directory('results/{sample}/salsa2'),
        fai = 'results/{sample}/salsa2/{sample}.fasta.fai',
        bed = 'results/{sample}/salsa2/{sample}.aln.bed'
    log:
        'logs/{sample}.salsa2_prep'
    conda:
        '../envs/salsa2_prep.yaml'
    threads:
        4
    resources:
        mem_mb = 48000
    params:
        runtime = '0-01:00:00'
    shell:
        "rm -rf {output.dir}; mkdir {output.dir};"
        
        "cp {input.assembly} {output.dir};"
        "bamToBed -i {input.aln} > {output.dir}/{wildcards.sample}.aln.bed;"
        
        "cd {output.dir};"
        "samtools faidx {wildcards.sample}.fasta;"
        "sort -k 4  {wildcards.sample}.aln.bed > tmp && mv tmp {wildcards.sample}.aln.bed;"
        
################################################################################
# salsa2                                                    #
################################################################################
        
rule salsa2:
    input:
        dir = 'results/{sample}/salsa2'
    output:
        'results/{sample}/salsa2.fa'
    log:
        'logs/{sample}.salsa2'
    conda:
        '../envs/salsa2.yaml'
    threads:
        8
    resources:
        mem_mb = 16000
    params:
        enzyme = lambda wildcards: samples_info.enzyme[wildcards.sample],
        runtime = '0-23:59:00'
    shell:
        "cp resources/run_pipeline.py $CONDA_PREFIX/bin/run_pipeline_{wildcards.sample}.py;"
        
        "chmod +x $CONDA_PREFIX/bin/run_pipeline_{wildcards.sample}.py;"
        
        "cd {input.dir}; rm -rf scaffolds;"

        "python $CONDA_PREFIX/bin/run_pipeline_{wildcards.sample}.py -a {wildcards.sample}.fasta \
        -l {wildcards.sample}.fasta.fai \
        -b {wildcards.sample}.aln.bed \
        -e {params.enzyme}  \
        -i 30 \
        -o scaffolds \
        -p yes \
        -m yes;"
        
        "cp scaffolds/scaffolds_FINAL.fasta ../salsa2.fa"
        
################################################################################
# salsa2_to_hic                                                    #
################################################################################
          
rule salsa2_to_hic:
    input:
        'results/{sample}/salsa2.fa'
    output:
        'results/{sample}/salsa2.hic'
    log:
        'logs/{sample}.salsa2_to_hic'
    conda:
        '../envs/convert.yaml'
    threads:
        4
    resources:
        mem_mb = 32000
    params:
        runtime = '0-12:00:00'
    shell:
        "bash resources/convert.sh \
        results/{wildcards.sample}/salsa2/scaffolds/;"
        "cp results/{wildcards.sample}/salsa2/scaffolds/salsa_scaffolds.hic {output};"
