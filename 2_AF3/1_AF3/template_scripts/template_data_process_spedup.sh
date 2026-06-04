#!/bin/bash
#SBATCH --job-name=mmcif_AF3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=genoa
#SBATCH --cpus-per-task=24
#SBATCH --time=01:00:00

#load Alphafold
module load 2024
module load AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0
# Avoids warnings
unset LD_PRELOAD

cd ${TMPDIR}
mkdir input_af3
mkdir output_af3

# Copy input jsons to tmpdir, direcotry input option
cp $1/*.json ./input_af3

# path to submission script AF3, single file option
#JSON_PATH=$1


#path to reduced database
DATA_PATH=/home/rvonk1/AF/modified_database

# path to modified script
python_file=/home/rvonk1/AF3_TCRpMHC_snellius/src/1_AF3/run_alphafold.py

# path to AF3 containing (fixed path for snellius surf)
AF3_CONTAINER_PATH=/sw/arch/RHEL9/EB_production/2024/software/AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0/bin/alphafold-3.0.0.sif

# path to directory of saved AF3 model weights
#model_root=/home/ddiepenbroek
#model_root=$2

# Set the proper project root
bind_root=/projects

#clock the time 
start_time=$(date +%s)
#
# single-file option
#cmd_args="--json_path $JSON_PATH 
# directory input option, running each json in the input directory
cmd_args="--input_dir=${TMPDIR}/input_af3
 --output_dir=${TMPDIR}/output_af3
 --db_dir ${DATA_PATH}
 --pdb_database_path ${DATA_PATH}/mmcif_files_AF3_selection
 --run_inference=False"
 # --model_dir ${model_root}

 # Probably does not work, run singularity exec instead.
 # alphafold-3.0.0.sif ${cmd_args}
 cmd="singularity exec --bind ${bind_root} ${AF3_CONTAINER_PATH} python ${python_file} ${cmd_args}"
 echo 'running command:'
 echo $cmd
 eval $cmd

# Copy output files back to home directory
cp -r ${TMPDIR}/output_af3/* $1/output_AF3_MSA/

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "The datapipeline was runned in $elapsed seconds"