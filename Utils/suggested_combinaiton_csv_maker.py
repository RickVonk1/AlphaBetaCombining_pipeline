#------------- Explanation -------------#
"""
Author : Rick Vonk

this makes a csv file with the full complex for the combinaitons that are a resulted promising combinaiton

"""
#------------ Import -------------#
import pandas as pd
from pathlib import Path
import csv
import os


#----------------- inputs ------------------#

ab_csv = '/home/rvonk1/3_Jolanda_data_Pipeline/Script_pipeline/1_pre-AF3/required_files/Alpha_Beta_single_chains.csv'

complex_list = ['a1_b12', 'a1_b8', 'a10_b8', 'a10_b11', 'a11_b12', 'a11_b4', 'a11_b6', 'a11_b10', 'a2_b7', 'a2_b12', 'a2_b11', 'a2_b1', 'a2_b2', 'a3_b11', 'a3_b8', 'a3_b7', 'a3_b1', 'a3_b12', 'a4_b8', 'a4_b7', 'a4_b12', 'a5_b10', 'a5_b4', 'a6_b8', 'a7_b12', 'a7_b8', 'a7_b10', 'a8_b11', 'a9_b11', 'a9_b8']

# ----------------- code ------------------ #


def df_maker(csv_file, complex_list):
    csv_df = pd.read_csv(csv_file)

    script_location = Path(__file__).resolve()
    script_dir = script_location.parent
    new_csv = os.path.join(script_dir,'csv_file_complex_sequence.csv')
    with open (new_csv, 'w', newline='') as f:
        writing = csv.writer(f)
        writing.writerow(['complex','tcr_a','tcr_b','peptide','mhc_a','mhc_b2m'])
        for entry in sorted(complex_list):
            tcr_a, tcr_b = entry.split('_')
            tcr_a_seq = csv_df.loc[csv_df['Number'] == tcr_a, 'Full sequence'].values[0]
            tcr_b_seq = csv_df.loc[csv_df['Number'] == tcr_b, 'Full sequence'].values[0]

            mhc_a_seq = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP'
            mhc_b_seq = 'IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM'
            peptide_seq = 'RLSSCVPV'
            writing.writerow([entry, tcr_a_seq, tcr_b_seq, peptide_seq, mhc_a_seq, mhc_b_seq])
        
    return new_csv

#----------------- Activation ------------------#
if __name__ == '__main__':
    df_maker(ab_csv, complex_list)
