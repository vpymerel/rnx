
################################################################################
# samtools sorting                                                          #
# necessary to use q option
################################################################################

rule samtools_sort:
    input:
        'results/{sample}/{sample}.bam'
    output:
        'results/{sample}/{sample}.sorted.bam'
    log:
        'logs/{sample}.samtools_sort'
    conda:
        '../envs/samtools.yaml'
    threads:
        8
    resources:
        mem_mb = 16000
    params:
        runtime = '0-01:00:00'
    shell:
        "samtools sort -n -@ 8 {input} -o {output};"

################################################################################
# yahs                                                          #
################################################################################
    
rule yahs:
    input:
        fasta = 'results/{sample}/mapped/{sample}.fasta',
        bam = 'results/{sample}/{sample}.sorted.bam'
    output:
        fasta = 'results/{sample}/yahs.fa',
        bin = 'results/{sample}/yahs/yahs.out.bin',
        agp = 'results/{sample}/yahs/yahs.out_scaffolds_final.agp',
    log:
        'logs/{sample}.yahs'
#    conda:
#        '../envs/yahs.yaml'
    threads:
        8
    resources:
        mem_mb = 400000
    params:
        enzyme = lambda wildcards: samples_info.enzyme[wildcards.sample],
        runtime = '0-03:00:00'
    shell:
        "mkdir -p results/{wildcards.sample}/yahs ;"
        
        "cd ./resources/yahs/; make -B; cd ../../;"
        "./resources/yahs/yahs -e {params.enzyme} \
        -q {config[QC]} \
        -r 500,750,1000,2500,5000,7500,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,20000000,50000000,100000000,200000000,500000000 \
        -o ./results/{wildcards.sample}/yahs/yahs.out \
        {input.fasta} {input.bam};"

        "cp results/{wildcards.sample}/yahs/yahs.out_scaffolds_final.fa {output.fasta}"

#
################################################################################
# juicer_pre                                                          #
################################################################################
            
rule juicer_pre:
    input:
        fai = 'results/{sample}/mapped/{sample}.fasta.fai',
        bam = 'results/{sample}/{sample}.sorted.bam',
        agp = 'results/{sample}/yahs/yahs.out_scaffolds_final.agp'
    output:
       txt = 'results/{sample}/yahs/out_JBAT.txt',
       log = 'results/{sample}/yahs/out_JBAT.log'
    log:
        'logs/{sample}.juicer_pre'
    conda:
        '../envs/yahs.yaml'
    threads:
        8
    resources:
        mem_mb = 32000
    params:
        runtime = '0-08:00:00',
        dir = 'results/{sample}/yahs/',
    shell:
        "juicer pre -a \
        -o {params.dir}/out_JBAT \
        {input.bam} \
        {input.agp} \
        {input.fai} >{output.log} 2>&1"

################################################################################
# yahs_to_hic                                                          #
################################################################################
 

rule yahs_to_hic:
    input:
       txt = 'results/{sample}/yahs/out_JBAT.txt',
       log = 'results/{sample}/yahs/out_JBAT.log'
    output:
       'results/{sample}/yahs.hic'
    log:
        'logs/{sample}.yahs_to_hic'
    conda:
        '../envs/juicer.yaml'
    threads:
        8
    resources:
        mem_mb = 32000
    params:
        runtime = '0-08:00:00',
        dir = 'results/{sample}/yahs/',
    shell:
        """
        (java -jar \
        -Xmx32G \
        resources/juicer_tools.1.9.9_jcuda.0.8.jar pre \
        {input.txt} \
        {output}.part <(cat {input.log}  | grep PRE_C_SIZE | awk '{{print $2" "$3}}')) && (mv {output}.part {output})
        """
     