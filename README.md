# Biobakery Data Processing Workflow 

A biobakery workflow pipeline for Whole Metagenome Shotgun (wmgx) data processing implemented using Snakemake rather than AnADAMA2. Programs supported thus far are Kneaddata, Metaphlan, and Humann2. 


![alt text](https://bitbucket.org/repo/5pd5AR/images/2528193080-wms_workflow.jpg)


https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home


## How to Run 
Warning - this is only to run what is in example-data. If real data is to be used Snakefile needs to be motified 

To run on HPC cluster:

    $ sh run_biobakery_workflow.sh
   
To run locally (for testing purposes) ##### *****WARNING computational heavy :

    $ /shares/hii/sw/snakemake/latest/bin/snakemake 



Singularity Image can be found at: 

    $ /shares/hii/images/morenoj/biobakery
