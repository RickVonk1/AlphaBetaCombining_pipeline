# !/bin/bash

python /home/rvonk1/3_Jolanda_data_Pipeline/Script_pipeline/2_AF3/2_postprocessing/tcr_af3_to_pdb_renumb.py \
    --input-dir /projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_3/2_AF/Process_2/TCR_post \
    --output-dir /projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_3/2_AF/Process_2/Renumbered_TCR

    # If you want to specify the path to the ImmunoPDB script of ANARCI, you can do so with the --immunopdb-path argument. Default is:
    # --immunopdb-path /projects/0/prjs1135/software/ANARCI/Example_scripts_and_sequences/ImmunoPDB.py