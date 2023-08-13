#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o ./out/%x.out
#SBATCH -e ./err/%x.err
##############################################
# USAGE: sbatch --job-name=prepForNetworkConstruction --cpus-per-task=16 --mem-per-cpu=4G prepForNetworkConstruction.sh $rds_object $name $out_dir
# Date 01/02/2023
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-R413

# Configure input arguments
rds_obj=$1
name=$2
out_dir=$3

echo -e "Loading rds object from: $rds_obj"
echo -e "Creating hdWGCNA in misc slot with name: $name"
echo -e "Outputting to: $out_dir\n"

# GRN inference
CMD="Rscript --vanilla prepForNetworkConstruction.R \
$rds_obj \
$name \
$out_dir"
echo -e "Running:\n $CMD\n"
$CMD

date