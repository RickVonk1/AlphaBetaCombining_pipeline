#!/bin/bash
#SBATCH --job-name=ensemble
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --partition=gpu_a100
#SBATCH --time=03:00:00
#SBATCH --error=./TCR_assembly_benchmark/4_log/ensemble/tcr_inf_a100_%j.err
#SBATCH --output=./TCR_assembly_benchmark/4_log/ensemble/tcr_inf_a100_%j.out

# load environment modules
module load 2024
module load AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0

cd /projects/0/prjs1135


#path to output directory
OUTPUT_PATH="./TCR_assembly_benchmark/3_docking_runs/output_AF3/test"


# Path of the container that is already hosted in the software stack of Snellius.
AF3_CONTAINER_PATH="/sw/arch/RHEL9/EB_production/2024/software/AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0/bin/alphafold-3.0.0.sif"

#path to the Alphafold3 weights
MODEL_PATH="/home/ddiepenbroek"

#clock the time 
start_time=$(date +%s)

# 
shopt -s globstar  # enable ** for recursive matching

for json_file in "$OUTPUT_PATH"/**/*_data.json; do
    echo "Processing: $json_file"

    # arguments
    cmd_args="--json_path ${json_file}
    --output_dir ${OUTPUT_PATH}
    --run_data_pipeline=False
    --model_dir ${MODEL_PATH}"

    unset LD_PRELOAD

    apptainer run --nv \
    -B "$PWD:/workspace" \
    --pwd /workspace \
    ${AF3_CONTAINER_PATH}  ${cmd_args}

done

#clock the end time
end_time=$(date +%s)
RUNTIME=$(( end_time - start_time ))
echo "Alphafold3 inference runtime: ${RUNTIME} seconds"
