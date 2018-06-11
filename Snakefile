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
METAPHLAN_PY = 'biobakery-metaphlan2-e7761e78f362/metaphlan2.py'
METAPHLAN_DIR ='biobakery-metaphlan2-e7761e78f362'

DEMO_SAMPLE = ['SRS014459-Stool.fasta.gz', 'SRS014464-Anterior_nares.fasta.gz', 'SRS014476-Supragingival_plaque.fasta.gz', 'SRS014494-Posterior_fornix.fasta.gz']

BOWTIE_GRCH38_BUILD = '/shares/hii/data/bt2/GRCh38'
NUCLEO_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/chocophlan'
PROTEIN_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/uniref'


rule all:
    input:
        expand('metaphlan_analysis/{samples}.txt', samples = DEMO_SAMPLE),
        expand('humann2_analysis/{samples}_genefamilies.tsv', samples = DEMO_SAMPLE)

#to be run only if preprocessing necssary (only accepts fastq files)
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
        expand('example-data/{samples}', samples =DEMO_SAMPLE)
    params:
        demo_data = DEMO_SAMPLE
    output:
        expand('metaphlan_analysis/{samples}.txt', samples = DEMO_SAMPLE)
    shell:
        """
        for sample in {params.demo_data}; do
            {SINGULARITY} python3 {METAPHLAN_PY} example-data/$sample \
                --input_type fasta \
                --bowtie2_exe {BOWTIE2_EXEC} \
                --bowtie2_build {BOWTIE2_BUILD} \
                --output metaphlan_analysis/$sample.txt
        done
        """

rule humann:
    input:
        expand('example-data/{samples}', samples =DEMO_SAMPLE)
    params:
        demo_data = DEMO_SAMPLE
    output:
        expand('humann2_analysis/{samples}_genefamilies.tsv', samples =DEMO_SAMPLE)
    shell:
        """
        for sample in {params.demo_data}; do
            {SINGULARITY} humann2 -i example-data/$sample \
                --bowtie2 {BOWTIE_DIR} \
                --metaphlan {METAPHLAN_DIR} \
                --taxonomic-profile metaphlan_analysis/$sample.txt \
                --nucleotide-database {NUCLEO_DB} \
                --protein-database {PROTEIN_DB} \
                -o humann2_analysis
        done
        """




