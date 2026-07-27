#------------- Explanation -------------#
"""
Author : Rick Vonk

This script makes a new experiment in the experiment folder and prepares required csv file for AF3 generation.
A target chain is chosen and combined with all beta chains

Input:
    - The Experiment directory
    - An experiment name
    - a alpha target chain from the csv file

Output:
    - Experiment folder
    - a csv file in the experiment folder with AF3 required data

Usage :
    python /path/to/script/
"""
#------------ Import -------------#
import csv
import shutil
import os
import re
import sys
from pathlib import Path

#----------------- inputs ------------------#
Exp_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

experiment = 'tcr_alpha_3'

target_a_chain = 'a3'

# ----------------- code ------------------ #
def Experiment_maker(output_dir, name):
    # Makes folders used for the AF3 pipeline
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir,'2_AF','Process_1'), exist_ok=True)
    os.makedirs(os.path.join(output_dir,'2_AF','Process_2','Raw_TCRs'), exist_ok=True)
    os.makedirs(os.path.join(output_dir,'2_AF','Process_2','Renumbered_TCR'), exist_ok=True)
    os.makedirs(os.path.join(output_dir,'2_AF','Process_2','TCR_post'), exist_ok=True)
    os.makedirs(os.path.join(output_dir,'Result'), exist_ok=True)
    os.makedirs(os.path.join(output_dir,'input_pdb'), exist_ok=True)
    return

def AF3_chains_maker(output_dir, chains_file, target):
    # Combines the TCR sequences wiht the given MHC sequence
    with open (chains_file,'r') as x:
        chains = csv.reader(x)
        header = next(chains)
        rows = list(chains)

    MHC = ["GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP","IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM","RLSSCVPV"]
    alpha_chain = []
    beta_chains = []
    for row in rows :
        if row[0] == target :
            alpha_chain.append(row)
        if 'b' in row[0] :
            beta_chains.append(row)

    path_target = os.path.join(output_dir,'2_AF/AF3_chains.csv')
    with open (path_target,'w',newline='') as x :
        AF3_csv = csv.writer(x)
        AF3_csv.writerow(["PDBID","A","B","M","N","P"])
        for case in beta_chains :
            AF3_csv.writerow([f"{target_a_chain}_{case[0]}",
                alpha_chain[0][11],
                case[11],
                MHC[0],
                MHC[1],
                MHC[2]
                ])
    return

#----------------- Activation ------------------#
working_dir = os.path.join(Exp_dir, experiment)

script_location = Path(__file__).resolve()
single_chains = os.path.join(script_location.parent, 'required_files', 'Alpha_Beta_single_chains.csv')


if __name__ == "__main__":
    Experiment_maker(working_dir, experiment)
    AF3_chains_maker(working_dir, single_chains, target_a_chain)



