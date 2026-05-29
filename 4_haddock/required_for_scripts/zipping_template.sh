#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=rome
#SBATCH --time=24:00:00
#SBATCH -o $Result_location/zipping/O.out
#SBATCH -e $Result_location/zipping/E.err


OUTPUT_FILE="$main_loc/$name.tar.gz"

tar -c --use-compress-program=pigz -f "$OUTPUT_FILE" "$haddock_output_loc" && rm -rf "$haddock_output_loc" && echo "Job Successful: $name has been zipped and removed."
