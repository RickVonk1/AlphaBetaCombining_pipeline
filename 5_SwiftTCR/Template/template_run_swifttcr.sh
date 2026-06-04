#!/bin/bash
#SBATCH --job-name=[case]_Swift
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=rome
#SBATCH --cpus-per-task=16
#SBATCH --time=04:30:00
#SBATCH --error=[output_dir]_log/run_%j.err
#SBATCH --output=[output_dir]_log/run_%j.out

START_TIME=$(date +%s)

# Activates conda enviroment
source activate swifttcr

# Goes into temporary space
cd ${TMPDIR}
mkdir [case]

echo "Path to output: ${TMPDIR}"

# Variable inputs, change to your own
INPUT_TCR_DIR="[tcr_loc]"
INPUT_MHC_DIR="[pmhc_loc]"
OUTPUT_DIR="[output_dir]"
TEMP_OUTPUT_DIR="${TMPDIR}/[case]"

SWIFTTCR_PATH="/home/rvonk1/swifttcr"
SWIFTTCR_PATH_PY="/home/rvonk1/swifttcr/scripts/swift_tcr.py"


# Docks every TCR in the input directory to every pMHC in the input directory
for mhc_file in "$INPUT_MHC_DIR"/*; do
    # Checks if the dir is not empty and contains a file
    echo "$mhc_file"
    [ -f "$mhc_file" ] || continue
    for tcr_file in "$INPUT_TCR_DIR"/*; do
        # Checks if the dir is not empty and contains a file
        echo "$tcr_file"
        [ -f "$tcr_file" ] || continue
        
        # Gets the tcr and cluster name
        base_tcr="${tcr_file##*/}"
        base_mhc="${mhc_file##*/}"
        mhc_name="${base_mhc%.*}" 
        name="${base_tcr%.*}_${mhc_name#*_}" 
        echo "$name"

        # Only runs swifttcr if tarred results do not exist yet
        if [ ! -f "${OUTPUT_DIR}/${name}/merged.tar" ]; then

            # Subshell to go to SwiftTCR and still being in TMPDIR
            ( cd ${SWIFTTCR_PATH} || { echo "Unable to change to Swifttcr dir"; exit 1; } 

            # Runs SwiftTCR
            python3 "$SWIFTTCR_PATH_PY" \
                -r "$mhc_file" \
                -l "$tcr_file" \
                -o "$TEMP_OUTPUT_DIR" \
                -op "$name" \
                -c 16 -t 3 -m 1000
            ) || { echo "Unable to run swifttcr for $name"; exit 1; } 

            # Tar merged dir and removes the dir
            tar -cvf "${TEMP_OUTPUT_DIR}/${name}/merged.tar" -C "${TEMP_OUTPUT_DIR}/${name}" merged && rm -r "${TEMP_OUTPUT_DIR}/${name}/merged" || { echo "tar failed for $name"; exit 1; } # && rm -r "$OUTPUT_DIR" <-- add after practice run and everything works fine
            mkdir -p "${OUTPUT_DIR}"

            # Copies the dir with merged.tar to given output dir
            cp -r "${TEMP_OUTPUT_DIR}/${name}" "${OUTPUT_DIR}" || { echo "Unable to copy the results for $name"; exit 1; } 
            echo "Directory is copied"
        else
            echo "${OUTPUT_DIR}/${name}/merged.tar exists"
        fi
    done
done


END_TIME=$(date +%s)
ELAPSED_SEC=$((END_TIME - START_TIME))

# Format seconds into HH:MM:SS
REMAINING_SEC=$((ELAPSED_SEC % 60))
TOTAL_MIN=$((ELAPSED_SEC / 60))
REMAINING_MIN=$((TOTAL_MIN % 60))
TOTAL_HOURS=$((TOTAL_MIN / 60))

printf "Total Execution Time: %02d:%02d:%02d\n" $TOTAL_HOURS $REMAINING_MIN $REMAINING_SEC