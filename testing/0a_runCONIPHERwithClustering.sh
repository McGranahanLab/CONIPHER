#!/bin/bash
#SBATCH --job-name=conipherRunVI             # Job name
#SBATCH --ntasks=1                         # Run a single task
#SBATCH --partition=cpu                    # Partition to run on
#SBATCH --mem=100G                          # Job Memory
#SBATCH --time=3-00:00:00                  # Time limit hrs:min:sec
#SBATCH --output=/nemo/project/proj-tracerx-wgs/working/tholk/pycloneVI_test/result_betabinomial/logs/conipherVIRun_%j.log   # Standard output and error log
#SBATCH --array=0-4                      # Range of indexes for array job and how many should run simultaneously

output_dir=/nemo/project/proj-tracerx-wgs/working/tholk/pycloneVI_test/result_betabinomial/

patients=(LTX0093 LTX0046 LTX0109 LTX0038)

mkdir -p ${output_dir}

sh /nemo/project/proj-tracerx-wgs/working/tholk/CONIPHER/testing/1a_runCONIPHERwithClustering.sh \
    ${patients[${SLURM_ARRAY_TASK_ID}]} \
    ${output_dir}
