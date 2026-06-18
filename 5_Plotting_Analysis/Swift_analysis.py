#------------- Explanation -------------#
"""
Author : Rick Vonk


"""
#------------ Import -------------#
import glob
import os
import pandas as pd

#----------------- inputs ------------------#
experiment_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

# ----------------- code ------------------ #

def data_extraction(experiment_dir):
    df_list = []
    search_path = os.path.join(experiment_dir, 'tcr_alpha_*')
    
    for chain in glob.glob(search_path):
        Swift_result = os.path.join(chain, 'Intermediate_files', 'Swift_Result.csv')
        
        df = pd.read_csv(Swift_result)
        df_list.append(df)
        
    df_combined = pd.concat(df_list, ignore_index=True)
    return df_combined

def csv_maker(experiment_dir, df):
    
    csv_loc = os.path.join(experiment_dir, 'Combined_results', 'Swift_Combined_data.csv')
    df.to_csv(csv_loc, index=False)

    return

#----------------- Activation ------------------#
if __name__ == '__main__':
    df = data_extraction(experiment_dir)
    csv_maker(experiment_dir, df)