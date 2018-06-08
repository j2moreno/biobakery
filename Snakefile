


import sys
import os

SINGULARITY_IMAGE = os.environ.get('SINGULARITY_IMAGE',
          '/shares/hii/images/morenoj/biobakery_image/test.simg')

SINGULARITY = '/shares/hii/sw/singularity/latest/bin/singularity exec '
SINGULARITY += '--bind /shares:/shares ' + SINGULARITY_IMAGE

BOWTIE2_EXEC = '/shares/hii/sw/bowtie2/2.2.9/bin/bowtie2'
BOWTIE2_BUILD = '/shares/hii/sw/bowtie2/2.2.9/bin/bowtie2-build'

BOWTIE_DIR = '/shares/hii/sw/bowtie2/2.2.9/bin'
TRIMMOMATIC_DIR = './bin'
METAPHLAN_DIR = './biobakery-metaphlan2-e7761e78f362/metaphlan2.py'

DEMO_SAMPLE = ['demo.fastq']
BOWTIE_GRCH38_BUILD = '/shares/hii/data/bt2/GRCh38'

rule all:
    input:
        'metaphlan_analysis/demo_metaphlan_output.txt'

rule kneaddata:
    input:
        expand('example-data/{samples}', samples = DEMO_SAMPLE)
    output:
        'test-kneaddata-output'

    shell:
        """
        {SINGULARITY} kneaddata \
            --input {input} \
            --reference-db {BOWTIE_GRCH38_BUILD} \
            --bowtie2 {BOWTIE_DIR} \
            --trimmomatic {TRIMMOMATIC_DIR} \
            --output {output}
        """


rule metaphlan:
    input:
        path = 'test-kneaddata-output',
        fastqFile = 'test-kneaddata-output/demo_kneaddata.fastq'
    output:
        'metaphlan_analysis/demo_metaphlan_output.txt'
    shell:
        """
        {SINGULARITY} python3 {METAPHLAN_DIR} {input.fastqFile} \
            --input_type fastq \
            --bowtie2_exe {BOWTIE2_EXEC} \
            --bowtie2_build {BOWTIE2_BUILD} \
            --output {output}
        """
