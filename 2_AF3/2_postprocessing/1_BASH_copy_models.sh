# !/bin/bash

# !!! Warning: this step will copy all the models from AF3 to a new directory. 
#      It should be changed in case you expect there to be too many models in the final directory.
INPUT_PATH=/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/Process_2/Raw_TCRs/test
OUTPUT_PATH=/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/Process_2/TCR_post

for case_path in "$INPUT_PATH"/*; do
    python /home/rvonk1/AF3_TCRpMHC_snellius/src/2_postprocessing/copy_and_name_models_af3.py --model-dir $case_path \
        --output-dir $OUTPUT_PATH \

    done
