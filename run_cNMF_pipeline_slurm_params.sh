#!/bin/bash

#SBATCH -J cnmf
#SBATCH -A pi-xuanyao
#SBATCH -p xuanyao-hm
#SBATCH --qos xuanyao
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=14-00:00:00
#SBATCH --error=./log_sbatch/cnmf_nf_dirver_%J.err
#SBATCH --output=./log_sbatch/cnmf_nf_dirver_%J.out

module load nextflow/24.10.6

# Limit NF driver
export NXF_OPTS='-Xms10G -Xmx20G'

# ################### #
#  invoke pipeline
# ################### #

nextflow run -resume main.nf \
    -profile hpc \
    -c nextflow_cnmf.config \
    -with-report reports_cnmf/reports_cnmf_$(date +%Y%m%d).html \
    -with-trace reports_cnmf/trace_cnmf_$(date +%Y%m%d).txt \
    -with-timeline reports_cnmf/timeline_$(date +%Y%m%d).html \
    -params-file cNMF_params.yml 

