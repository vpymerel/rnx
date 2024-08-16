################################################################################
# cfg_juicer                                                             #
################################################################################

rule cfg_juicer:
    input:
        assembly = 'data/{sample}.fasta',
        R1 = 'results/{sample}/{sample}_R1.fastq.gz',
        R2 = 'results/{sample}/{sample}_R2.fastq.gz'
    output:
        dir = directory('results/{sample}/juicer/'),
        fa = 'results/{sample}/juicer/references/{sample}.fasta'

    log:
        'logs/{sample}.cfg_juicer'
    conda:
        '../envs/bwa.yaml'
    threads:
        2
    resources:
        mem_mb = 8000
    params:
        runtime = '0-01:00:00'
    shell:
        "rm -rf {output.dir}; mkdir {output.dir};"
        "cp -r resources/juicer/* {output.dir};"
        "cd {output.dir};"
        "ln -s CPU scripts;"
        "cd scripts/common;"
        "wget https://github.com/aidenlab/Juicebox/releases/download/v2.17.00/juicer_tools_2.17.00.jar;"
        "ln -s juicer_tools_2.17.00.jar juicer_tools.jar;"
        "cd ../.. ;"
        "mkdir references;"
        "cd references;"
        "cp ../../../../{input.assembly} ./ ;"
        "bwa index *fasta;"
        "cd ../ ;"
        "mkdir work;"
        "cd work;"
        "mkdir fastq;"
        "cd fastq;"
        "cp ../../../../../{input.R1} ./{wildcards.sample}_R1.fastq.gz;"
        "cp ../../../../../{input.R2} ./{wildcards.sample}_R2.fastq.gz;"

################################################################################
# cs_size                                                             #
################################################################################

rule cs_size:
    input:
        'results/{sample}/juicer/references/{sample}.fasta'
    output:
         'results/{sample}/juicer/references/{sample}.sizes'
    log:
        'logs/{sample}.cs_size'
    conda:
        '../envs/bioawk.yaml'
    threads:
        2
    resources:
        mem_mb = 2000
    params:
        runtime = '0-00:10:00'
    shell:
        """
        bioawk -c fastx '{{print $name"\\t"length($seq)}}' {input} > {output};
        """

################################################################################
# restriction_sites                                                             #
################################################################################
        
rule restriction_sites:
    input:
        assembly = 'results/{sample}/juicer/references/{sample}.fasta'
    output:
         'results/{sample}/juicer/references/{sample}.restriction_sites'
    log:
        'logs/{sample}.restriction_sites'
    conda:
        '../envs/restriction_sites.yaml'
    threads:
        2
    resources:
        mem_mb = 2000
    params:
        enzyme = get_sample_enzyme_as_list,
        runtime = '0-03:00:00'
    shell:
        """
        sed -e "s/'ToRemplace'/['GATC', 'GANTC', 'CTNAG', 'TTAA']/g" resources/juicer/misc/generate_site_positions.py> resources/juicer/misc/generate_site_positions.{wildcards.sample}.py;
        """
        "python3 resources/juicer/misc/generate_site_positions.{wildcards.sample}.py \
        UserDefined \
        {wildcards.sample} \
        {input.assembly};"
        "mv {wildcards.sample}_UserDefined.txt {output}"



################################################################################
# juicer                                                             #
################################################################################
  
rule juicer:
    input:
        dir = 'results/{sample}/juicer/',
        assembly = 'results/{sample}/juicer/references/{sample}.fasta',
        assembly_sizes = 'results/{sample}/juicer/references/{sample}.sizes',
        restriction_sites = 'results/{sample}/juicer/references/{sample}.restriction_sites'
    output:
        txt = 'results/{sample}/juicer/work/aligned/merged_nodups.txt',
    log:
        'logs/{sample}.juicer'
    conda:
        '../envs/juicer.yaml'
    threads:
        24
    resources:
        mem_mb = 48000
    params:
        runtime = '0-24:00:00'
    shell:
        "rm -rf {input.dir}/work/splits;"
        "rm -rf {input.dir}/work/aligned;"
        "{input.dir}/scripts/juicer.sh \
        -g {wildcards.sample} \
        -d {config[folder]}/{input.dir}/work \
        -a {wildcards.sample} \
        -p {config[folder]}/{input.assembly_sizes} \
        -y {config[folder]}/{input.restriction_sites} \
        -z {config[folder]}/{input.assembly} \
        -D {config[folder]}/{input.dir} \
        -t {threads} \
        -S early \
        -f;"   


#        -s {params.enzyme} \
