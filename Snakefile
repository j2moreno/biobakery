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

DEMO_SAMPLES = ['131014_SN805_0422_AC2UG1ACXX_5_IDMB9_', '131014_SN805_0422_AC2UG1ACXX_5_IDMB92_', '130525_SN7001324_0153_AD253AACXX_2_IDMB81_', \
        '130525_SN7001324_0153_AD253AACXX_7_IDMB48_', '130525_SN7001324_0153_AD253AACXX_7_IDMB50_', '130525_SN7001324_0153_AD253AACXX_7_IDMB68_']
KNEADDATA_LOGS = 'output/kneaddata_output/logs'

BOWTIE_GRCH38_BUILD = '/shares/hii/data/bt2/GRCh38'
NUCLEO_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/chocophlan'
PROTEIN_DB = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/uniref'
UTILITY_MAPPING = '/shares/hii/bioinfo/ref/biobakery-workflows/humann2/utility_mapping/map_level4ec_uniref90.txt.gz'

rule all:
    input:
        expand('output/kneaddata_output/{sample}/{sample}R1_R2_.fastq', sample = DEMO_SAMPLES),
        'vis_inputs/kneaddata/merged/kneaddata_read_counts.tsv',
        expand('output/metaphlan_analysis/{sample}taxo_bugs.txt', sample = DEMO_SAMPLES),
        'vis_inputs/metaphlan2/merged/metaphlan2_taxonomic_profiles.tsv',
        #expand('output/metaphlan_analysis/{sample}*.sam', sample = DEMO_SAMPLES),
        expand('output/humann2_analysis/{sample}/{sample}R1_R2__genefamilies.tsv', sample = DEMO_SAMPLES),
        'vis_inputs/humann2/counts/humann2_read_and_species_count_table.tsv'

#to be run only if preprocessing necssary (only accepts fastq files)
rule kneaddata:
    input:
        R1 = 'example-data/{sample}1.fastq',
        R2 = 'example-data/{sample}2.fastq'
    params:
        directory = 'output/kneaddata_output/{sample}'
    output:
        'output/kneaddata_output/{sample}/{sample}R1_R2_.fastq'
    shell:
        """
        {SINGULARITY} kneaddata \
            --input {input.R1} \
            --input {input.R2} \
            --reference-db {BOWTIE_GRCH38_BUILD} \
            --bowtie2 {BOWTIE_DIR} \
            --trimmomatic {TRIMMOMATIC_DIR} \
            --output {params.directory}

        cat {params}/*_paired_*.fastq >> {output}

        """

rule kneaddata_post_processing:
    input:
        expand('output/kneaddata_output/{sample}/{sample}R1_R2_.fastq', sample = DEMO_SAMPLES)
    output:
        'vis_inputs/kneaddata/merged/kneaddata_read_counts.tsv'
    shell:
        """
        find output/kneaddata_output/ -name "*.log" | while read i; do \
            mv $i output/kneaddata_output/logs/; done

        {SINGULARITY} kneaddata_read_count_table \
            --input {KNEADDATA_LOGS} \
            --output {output}
        """


rule metaphlan:
    input:
        'output/kneaddata_output/{sample}/{sample}R1_R2_.fastq'
    output:
        'output/metaphlan_analysis/{sample}taxo_bugs.txt'
    shell:
        """
        {SINGULARITY} python {METAPHLAN_PY} {input} \
            --input_type fastq \
            --bowtie2_exe {BOWTIE2_EXEC} \
            --bowtie2_build {BOWTIE2_BUILD} \
            --samout output/metaphlan_analysis/$(basename {input}).sam \
            --output {output}
        """

rule metaphlan_post_processing:
    input:
        expand('output/metaphlan_analysis/{sample}taxo_bugs.txt', sample = DEMO_SAMPLES)
    output:
        'vis_inputs/metaphlan2/merged/metaphlan2_taxonomic_profiles.tsv'

    shell:
        """
        {SINGULARITY} python {METAPHLAN_DIR}/utils/merge_metaphlan_tables.py {input} \
            > {output}
        """

rule humann:
    input:
        fastq = 'output/kneaddata_output/{sample}/{sample}R1_R2_.fastq',
        profile = 'output/metaphlan_analysis/{sample}taxo_bugs.txt',
    params:
        directory = 'output/humann2_analysis/{sample}'
    output:
        'output/humann2_analysis/{sample}/{sample}R1_R2__genefamilies.tsv'
    shell:
        """
        {SINGULARITY} humann2 -i {input.fastq} \
            --bowtie2 {BOWTIE_DIR} \
            --metaphlan {METAPHLAN_DIR} \
            --taxonomic-profile {input.profile} \
            --nucleotide-database {NUCLEO_DB} \
            --protein-database {PROTEIN_DB} \
            --threads 4 \
            -o {params}
        """


rule huamnn_post_processing:
    input:
        expand('output/humann2_analysis/{sample}/{sample}R1_R2__genefamilies.tsv', sample = DEMO_SAMPLES)
    params:
        'output/humann2_analysis/path_Abundance_Norm'
    output:
        'vis_inputs/humann2/counts/humann2_read_and_species_count_table.tsv'
    shell:
        """
        find output/humann2_analysis/ -name "*_pathabundance.tsv" | while read i; do \
            echo $i
            echo $(basename $i)
            {SINGULARITY} humann2_renorm_table -i $i \
                --units relab \
                --output {params}/$(basename $i); \
        done

        {SINGULARITY} humann2_join_tables \
            --input output/humann2_analysis/path_Abundance_Norm \
            --output output/humann2_analysis/path_Abundance_Norm/pathAbun_merged_tables.tsv


        find output/humann2_analysis/ -name "*_genefamilies.tsv" | while read i; do \
            {SINGULARITY} humann2_regroup_table -i $i \
                -c {UTILITY_MAPPING} \
                --output output/humann2_analysis/geneFamilies_Grouped/$(basename $i)_ecs.tsv; \
        done

        find output/humann2_analysis/ -name "*_genefamilies.tsv" | while read i; do \
            {SINGULARITY} humann2_renorm_table -i $i \
                --units relab \
                --output output/humann2_analysis/geneFamilies_Norm/$(basename $i)_genefamilies_norm.tsv; \
        done



        find output/humann2_analysis/ -name "*ecs.tsv" | while read i; do \
            {SINGULARITY} humann2_renorm_table -i $i \
                --units relab \
                --output output/humann2_analysis/geneFamilies_Grouped_Norm/$(basename $i)_ecs_norm.tsv; \
        done

        {SINGULARITY} humann2_join_tables \
            --input output/humann2_analysis/geneFamilies_Grouped_Norm \
            --output output/humann2_analysis/geneFamilies_Grouped_Norm/EC_merged.tsv



        {SINGULARITY} humann2_join_tables \
            --input output/humann2_analysis/geneFamilies_Norm \
            --output output/humann2_analysis/geneFamilies_Norm/geneFam_merged_tables.tsv


        {SINGULARITY} python biobakery-scripts/count_features.py \
            -i output/humann2_analysis/geneFamilies_Norm/geneFam_merged_tables.tsv \
            --reduce-sample-name --ignore-un-features --ignore-stratification \
            -o vis_inputs/humann2/merged/geneFamilies_counts.tsv

        {SINGULARITY} python biobakery-scripts/count_features.py \
            -i output/humann2_analysis/geneFamilies_Grouped_Norm/EC_merged.tsv \
            --reduce-sample-name --ignore-un-features --ignore-stratification \
            -o vis_inputs/humann2/merged/EC_counts.tsv

        {SINGULARITY} python biobakery-scripts/count_features.py \
            -i output/humann2_analysis/path_Abundance_Norm/pathAbun_merged_tables.tsv \
            --reduce-sample-name --ignore-un-features --ignore-stratification \
            -o vis_inputs/humann2/merged/pathabundance_relab.tsv

        {SINGULARITY} humann2_join_tables -i vis_inputs/humann2/merged \
            --file_name _counts.tsv \
            -o vis_inputs/humann2/counts/humann2_feature_counts.tsv

        find output/humann2_analysis/ -name "*_pathcoverage.tsv" | while read i; do \
            cp $i output/humann2_analysis/path_Coverage/; done

        touch vis_inputs/anadama.log

        find output/humann2_analysis/ -name "*.log" | while read i; do \
            cp $i output/humann2_analysis/logs/; done

        {SINGULARITY} python biobakery-scripts/get_counts_from_humann2_logs.py \
            --input output/humann2_analysis/logs \
            --output vis_inputs/humann2/counts/humann2_read_and_species_count_table.tsv

        """
