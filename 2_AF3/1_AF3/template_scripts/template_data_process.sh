#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=genoa
#SBATCH --cpus-per-task=24
#SBATCH --time=13:00:00

module load 2024
module load AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0

cd ${TMPDIR}
mkdir input_af3
mkdir output_af3

echo "Path to output: ${TMPDIR}"
echo ${ls} 

# Copy input jsons to tmpdir
cp $1/*.json ./input_af3

# Path of the container that is already hosted in the software stack of Snellius.
AF3_CONTAINER_PATH="/sw/arch/RHEL9/EB_production/2024/software/AlphaFold/3.0.0-foss-2024a-CUDA-12.6.0/bin/alphafold-3.0.0.sif"

# Path of the data. Contains both the (large) .fasta files used for jackhmmer and a simlink to the mmcif files on NVME storage used for MSA deduplication and template matching.
DATA_PATH=/projects/2/managed_datasets/AlphaFold/3.0.0

# Path to the mmcif symlink. For now we point directly to the 'true' location at /scratch-nvme/ml-datasets/AlphaFold/3.0.0/mmcif_files/.
# Setting it to {DATA_PATH}/mmcif_files should also work due to the symlink
MMCIF_PATH=/scratch-nvme/ml-datasets/AlphaFold/3.0.0/mmcif_files/

# AF3 command line arguments
cmd_args="--input_dir=${TMPDIR}/input_af3
--output_dir=${TMPDIR}/output_af3
--db_dir=${DATA_PATH}
--pdb_database_path=${MMCIF_PATH}
--run_inference=False" # Do not run inference

# Unset to avoid warnings.
unset LD_PRELOAD

# Record start time
START=$(date +%s)

# Run the Alphafold 3 data pipeline.
# -B "$PWD:/workspace" mounts the current directory ($PWD) to /workspace inside the container.
# -B ${DATA_PATH} mounts the data path to the container.
# --pwd sets the working directory inside the container.
apptainer run -B "$PWD:/workspace" \
    -B ${DATA_PATH} \
    --pwd /workspace \
    ${AF3_CONTAINER_PATH}  ${cmd_args}

# Record end time
END=$(date +%s)

# Compute runtime
RUNTIME=$((END - START))
echo "Alphafold3 runtime: ${RUNTIME} seconds"

# Copy output files back to home directory
cp -r ${TMPDIR}/output_af3/* $1/output_AF3_MSA/