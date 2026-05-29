# Activate the conda environment if needed
# source activate alphafold3

# Submit AF3 inference jobs
python /home/rvonk1/AF3_TCRpMHC_snellius/src/1_AF3/preprocess_and_run.py --mode run_inference \
    --input-csv /projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/AF3_chains.csv \
    --output-dir /projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/tcr_alpha_11/2_AF/Process_1 \
    --model-weights /home/rvonk1/AF/Weights \
    --tcrpmhc-specific \
    --chainID-columns A B M N P #Column ID has to be pdb-standard, one-letter uppercase