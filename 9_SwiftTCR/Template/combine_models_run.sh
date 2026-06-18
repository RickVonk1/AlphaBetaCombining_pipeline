#!/bin/bash
#SBATCH --job-name=[case]_combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=genoa
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH --output=[output_dir]_log/Combining_%j.out
#SBATCH --error=/[output_dir]_log/Combining_%j.err

set -euo pipefail

source activate swifttcr
cd /home/rvonk1/swifttcr

# Zorg dat de python output direct geprint wordt naar je log
export PYTHONUNBUFFERED=1

base_dir="[output_dir]"

echo "Start clustering"
complete_start=$(date +%s)

# for dir in "$base_dir"/*; do
#   # Record start time
#   start_pdb=$(date +%s)
#   echo "$dir"
python3 [python_script] \
  --input_dir [output_dir] \
  --output_dir [combined_dir] \
  --cores ${SLURM_CPUS_PER_TASK}
  # Record end time
  # end_pdb=$(date +%s)

#   # Compute runtime
#   runtime_pdb=$((end_pdb - start_pdb))
#   echo "PDB runtime: ${runtime_pdb}"
# done

complete_end=$(date +%s)
complete_runtime=$((complete_end - complete_start))
echo "Runtime pdbName: ${complete_runtime}"
echo "End of job"

