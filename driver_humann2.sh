#!/usr/bin/env bash

sbatch \
    --partition=hii02 \
    --cpus-per-task=4 \
    --mem=20G \
    --time=1-0 \
    --output=test.log \
    run_Humann2.sh \

