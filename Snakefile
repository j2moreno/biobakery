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

DEMO_SAMPLES = ['SRS014459-Stool', 'SRS014464-Anterior_nares', 'SRS014476-Supragingival_plaque', 'SRS014494-Posterior_fornix']

BOWTIE_GRCH38_BUILD = '/shares/hii/data/bt2/GRCh38'
NUCLEO_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/chocophlan'
PROTEIN_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/uniref'


rule all:
    input:
        expand('output/metaphlan_analysis/{sample}.txt', sample = DEMO_SAMPLES),
        #expand('humann2_analysis/{sample}_genefamilies.tsv', sample = DEMO_SAMPLES)

#to be run only if preprocessing necssary (only accepts fastq files)
rule kneaddata:
    input:
        'example-data/{sample}.fasta.gz'
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
        'example-data/{sample}.fasta.gz'
    output:
        'output/metaphlan_analysis/{sample}.txt'
    shell:
        """
        {SINGULARITY} python3 {METAPHLAN_PY} {input} \
            --input_type fasta \
            --bowtie2_exe {BOWTIE2_EXEC} \
            --bowtie2_build {BOWTIE2_BUILD} \
            --bowtie2out {output}.bowtie2out.txt \
            --output {output}
        """

#rule humann:
#    input:
#        expand('metaphlan_analysis/{samples}.txt', samples = DEMO_SAMPLE)
#    params:
#        demo_data = DEMO_SAMPLE
#    output:
#        expand('humann2_analysis/{samples}_genefamilies.tsv', samples =DEMO_SIMPLE)
#    shell:
#        """
#        for sample in {params.demo_data}; do
#            {SINGULARITY} humann2 -i example-data/$sample \
#                --bowtie2 {BOWTIE_DIR} \
#                --metaphlan {METAPHLAN_DIR} \
#                --taxonomic-profile metaphlan_analysis/$sample.txt \
#                --nucleotide-database {NUCLEO_DB} \
#                --protein-database {PROTEIN_DB} \
#                -o humann2_analysis
#        done
#        """
#



