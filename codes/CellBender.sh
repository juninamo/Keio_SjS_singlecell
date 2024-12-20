#!/bin/sh
#$ -S /bin/sh
#SBATCH --array=1-27
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4

# conda create -n cellbender python=3.7
module load anaconda
conda activate cellbender
# pip install cellbender

seq_libs=(List LB094 LB098 LB099 LB100 LB101 LB103 LB104 LB105 LB106 LB107 LB110 LB111 LB116 LB118 LB121 LB125 LB127 LB128 LB130 LB132 LB133 LB134 LB135 LB138 LB139 LB140 LB142)
sample_id=${seq_libs[$SLURM_ARRAY_TASK_ID]}

DIR=/path/to/
cd ${DIR}/${sample_id}/

inputh5="raw_feature_bc_matrix.h5"
outputh5="cellbender_feature_bc_matrix.h5"

cellbender remove-background \
        --input ${inputh5} \
        --output ${outputh5} \
        --fpr 0.01 \
        --epochs 150 \
        --cpu-threads 4 \
        --debug
