################################################################################
# ThreeD_DNA                                                    #
################################################################################

rule ThreeD_DNA:
    input:
        assembly = 'results/{sample}/juicer/references/{sample}.fasta',
        juicer = 'results/{sample}/juicer/work/aligned/merged_nodups.txt'
    output:
        dir = directory('results/{sample}/ThreeD_DNA'),
        fasta = 'results/{sample}/ThreeD_DNA.fa',
        hic = 'results/{sample}/ThreeD_DNA.hic'
    log:
        'logs/{sample}.ThreeD_DNA'
    conda:
        '../envs/ThreeD_DNA.yaml'
    threads:
        8
    resources:
        mem_mb = 80000
    params:
        runtime = '0-23:59:00'
    shell:
        "rm -rf {output.dir}; mkdir {output.dir};"
        "cp -r resources/3D_DNA/* {output.dir};"
        "cp {input.assembly} {output.dir};"
        "cp {input.juicer} {output.dir};"
        "cd {output.dir};"
        "bash run-asm-pipeline.sh -q {config[QC]} {wildcards.sample}.fasta merged_nodups.txt;"
        "cp {wildcards.sample}.final.fasta ../../../{output.fasta};"
        "cp {wildcards.sample}.final.hic ../../../{output.hic};"
