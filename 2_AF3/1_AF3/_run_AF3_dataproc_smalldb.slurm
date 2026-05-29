#!/bin/bash
#SBATCH --job-name=J_data_AF3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=genoa
#SBATCH --cpus-per-task=24
#SBATCH --time=01:00:00
#SBATCH --error=/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/Process_2/log/data/%j.err
#SBATCH --output=/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/Process_2/log/data/%j.out

#load Alphafold
module load 2024
module load AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0
# Avoids warnings
unset LD_PRELOAD

#path to reduced database
data_root=/home/rvonk1/AF/modified_database

# path to modified script
python_file=/home/rvonk1/3_Jolanda_data_Pipeline/Script_pipeline/2_AF3/1_AF3/run_alphafold.py

# path to AF3 containing (fixed path for snellius surf)
container_root=/sw/arch/RHEL9/EB_production/2024/software/AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0/bin/alphafold-3.0.0.sif

# path to directory of saved AF3 model weights
model_root=/home/rvonk1/AF/Weights

# path to submission script AF3
JSON_PATH="$1"

# Set the proper project root
project_root=/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/Process_2/Raw_TCRs
bind_root=/projects

#clock the time 
start_time=$(date +%s)
#
cmd_args="--json_path $JSON_PATH
 --output_dir ${project_root}/test/
 --db_dir ${data_root}
 --pdb_database_path ${data_root}/mmcif_files_AF3_selection
 --model_dir ${model_root}
 --run_inference=False"

 # Probably does not work, run singularity exec instead.
 # alphafold-3.0.0.sif ${cmd_args}
 cmd="singularity exec --bind ${bind_root} ${container_root} python ${python_file} ${cmd_args}"
 echo 'running command:'
 echo $cmd
 eval $cmd

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "The datapipeline was runned in $elapsed seconds"