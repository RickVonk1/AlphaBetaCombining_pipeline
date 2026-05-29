#!/bin/bash
#SBATCH --job-name=$iteration_haddock3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=genoa
#SBATCH --cpus-per-task=24
#SBATCH --time=10:00:00
#SBATCH --output=$working_dir/$exp_name/Haddock_output/log/_%j.out
#SBATCH --error=$working_dir/$exp_name/Haddock_output/log/_%j.err

source activate haddock3


$cases_to_run

wait
