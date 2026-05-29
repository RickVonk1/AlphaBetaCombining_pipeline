#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=rome
#SBATCH --time=24:00:00
#SBATCH -o $Result_location/zipping/O.out
#SBATCH -e $Result_location/zipping/E.err


OUTPUT_FILE="$main_loc/$name.tar.gz"

cd "$main_loc" || exit 1

tar -c --use-compress-program=pigz -f "$OUTPUT_FILE" "$folder_to_zip" && rm -rf "$folder_to_zip"