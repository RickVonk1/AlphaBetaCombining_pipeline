#!/bin/bash

# This will given a list of json files loop over said json files and run the two AF3 scripts


JSON_DIR=/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/Process_1

for JSON_PATH in "$JSON_DIR"/*/*.json
do

    base_name=$(basename "$JSON_PATH" .json)

    pdb_part=${base_name:0:4}
    rest_part=${base_name:4}
    base_name="$(echo "$pdb_part" | tr '[:lower:]' '[:upper:]')$rest_part"

    echo "Submitting pipeline for $base_name"

    # Submit first job and capture job ID
    job1_id=$(sbatch \
        --job-name=AF3_${base_name} \
        /home/rvonk1/3_Jolanda_data_Pipeline/Script_pipeline/2_AF3/1_AF3/_run_AF3_dataproc_smalldb.slurm "$JSON_PATH" \
        | awk '{print $4}')

    # Submit second job dependent on first finishing successfully
    sbatch \
        --job-name=Inf_${base_name} \
        --dependency=afterok:$job1_id \
        /home/rvonk1/3_Jolanda_data_Pipeline/Script_pipeline/2_AF3/1_AF3/_run_AF3_inf_smalldb.slurm "$base_name"

done