#!/bin/bash

patient=$1
shift
output_dir=$1

# cd /nemo/project/proj-tracerx-lung/tctProjects/mcgranahanLab/tholk/apps/pyclone-vi/
# mamba create -n conipherVI -c conda-forge -c bioconda conipher --file requirements.txt --yes

source activate /camp/home/tholk/.conda/envs/conipherVI

Rscript /nemo/project/proj-tracerx-wgs/working/tholk/CONIPHER/testing/2a_runCONIPHERwithClustering.R ${patient} \
    ${output_dir}/${patient}/

conda deactivate