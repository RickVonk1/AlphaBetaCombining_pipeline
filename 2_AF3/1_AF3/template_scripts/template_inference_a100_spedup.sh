#!/bin/bash
#SBATCH --job-name=mmcif_AF3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gpu_a100
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=18
#SBATCH --time=01:00:00
#SBATCH --error=/projects/0/prjs1135/report_danielle/AF3_modified_database/test/small_database_test_inf.err
#SBATCH --output=/projects/0/prjs1135/report_danielle/AF3_modified_database/test/small_database_test_inf.out


module load 2024
module load AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0
# Avoids warnings
unset LD_PRELOAD

#path to reduced database
DATA_PATH=/home/rvonk1/AF/modified_database

# path to modified script
python_file=/home/rvonk1/AF3_TCRpMHC_snellius/src/1_AF3/run_alphafold.py

# path to AF3 containing (fixed path for snellius surf)
AF3_CONTAINER_PATH=/sw/arch/RHEL9/EB_production/2024/software/AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0/bin/alphafold-3.0.0.sif

# path to directory of saved AF3 model weights
#model_root=/home/ddiepenbroek
model_root=$1

# path to submission script AF3
INPUT_PATH=$2
#/projects/0/prjs1135/report_danielle/AF3_modified_database/test/8shi_tcrmodel2_data/8shi_tcrmodel2_data_data.json

# Set the proper project root
bind_root=/projects

#clock the time 
start_time=$(date +%s)
#
for JSON_PATH in ${INPUT_PATH}/output_AF3_MSA/*_rs*/*.json; do
    echo "Processing: $JSON_PATH"
    cmd_args=" 
    --json_path $JSON_PATH
    --output_dir ${INPUT_PATH}/AF3_inference_output/
    --db_dir ${DATA_PATH}
    --pdb_database_path ${DATA_PATH}/mmcif_files_AF3_selection
    --run_data_pipeline=False
    --model_dir ${model_root}"
    
    # Probably does not work, run singularity exec instead.
    # alphafold-3.0.0.sif ${cmd_args}
    cmd="singularity exec --nv --bind ${bind_root} ${AF3_CONTAINER_PATH} python ${python_file} ${cmd_args}"
    echo 'running command:'
    echo $cmd
    eval $cmd
done

#clock the end time
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "The inference pipeline was runned in $elapsed seconds"