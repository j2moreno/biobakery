import sys
import os

SINGULARITY_IMAGE = os.environ.get('SINGULARITY_IMAGE',
          '/shares/hii/images/morenoj/biobakery/test.simg')

SINGULARITY = '/shares/hii/sw/singularity/latest/bin/singularity exec '
SINGULARITY += '--bind /shares:/shares ' + SINGULARITY_IMAGE

BOWTIE2_EXEC = '/shares/hii/sw/bowtie2/2.2.9/bin/bowtie2'
BOWTIE2_BUILD = '/shares/hii/sw/bowtie2/2.2.9/bin/bowtie2-build'

BOWTIE_DIR = '/shares/hii/sw/bowtie2/2.2.9/bin'
TRIMMOMATIC_DIR = './bin'
METAPHLAN_PY = 'biobakery-metaphlan2-e7761e78f362/metaphlan2.py'
METAPHLAN_DIR ='biobakery-metaphlan2-e7761e78f362'

DEMO_SAMPLES = ['131014_SN805_0422_AC2UG1ACXX_5_IDMB9_', '131014_SN805_0422_AC2UG1ACXX_5_IDMB92_']
KNEADDATA_LOGS = 'output/kneaddata_output/logs'

BOWTIE_GRCH38_BUILD = '/shares/hii/data/bt2/GRCh38'
NUCLEO_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/chocophlan'
PROTEIN_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/uniref'


rule all:
    input:
        expand('output/kneaddata_output/{sample}', sample = DEMO_SAMPLES),
        'output/kneaddata_output/logs/kneaddata_read_counts.txt',
        expand('output/metaphlan_analysis/{sample}/{sample}taxo_bugs.txt', sample = DEMO_SAMPLES),
        #expand('humann2_analysis/{sample}_genefamilies.tsv', sample = DEMO_SAMPLES)

#to be run only if preprocessing necssary (only accepts fastq files)
rule kneaddata:
    input:
        R1 = 'example-data/{sample}1.fastq',
        R2 = 'example-data/{sample}2.fastq'
    output:
        'output/kneaddata_output/{sample}/merged_R1_R2_.fastq'
    shell:
        """
        {SINGULARITY} kneaddata \
            --input {input.R1} \
            --input {input.R2} \
            --reference-db {BOWTIE_GRCH38_BUILD} \
            --bowtie2 {BOWTIE_DIR} \
            --trimmomatic {TRIMMOMATIC_DIR} \
            --output {output}

        cat {output}/*_paired_*.fastq >> {output}

        """
rule kneaddata_post_processing:
    input:
        expand('output/kneaddata_output/{sample}/merged_R1_R2_.fastq', sample = DEMO_SAMPLES)
    output:
        'output/kneaddata_output/logs/kneaddata_read_counts.txt'
    shell:
        """
        find output/kneaddata_output/ -name "*.log" | while read i; do \
            cp $i output/kneaddata_output/logs/; done

        {SINGULARITY} kneaddata_read_count_table \
            --input {KNEADDATA_LOGS} \
            --output {output}
        """


rule metaphlan:
    input:
        'output/kneaddata_output/{sample}/merged_R1_R2_.fastq'
    output:
        'output/metaphlan_analysis/{sample}/{sample}taxo_bugs.txt'
    shell:
        """
        {SINGULARITY} python {METAPHLAN_PY} {input} \
            --input_type fastq \
            --bowtie2_exe {BOWTIE2_EXEC} \
            --bowtie2_build {BOWTIE2_BUILD} \
            --bowtie2out output/metaphlan_analysis/{sample}/{sample}.bowtie2out.txt \
            --output {output}
        """

rule humann:
    input:
        'output/metaphlan_analysis/{sample}.txt'
    output:
        'output/humann2_analysis/{sample}_genefamilies.tsv'
    shell:
        """
        {SINGULARITY} humann2 -i {input} \
            --bowtie2 {BOWTIE_DIR} \
            --metaphlan {METAPHLAN_DIR} \
            --taxonomic-profile output/metaphlan_analysis/{sample}.txt \
            --nucleotide-database {NUCLEO_DB} \
            --protein-database {PROTEIN_DB} \
            -o humann2_analysis
        """




