#!/usr/bin/env bash

sbatch \
    --partition=hii02 \
    --cpus-per-task=4 \
    --mem=60G \
    --time=1-0 \
    --output=test.log \
    run_biobakery_workflow.sh \

