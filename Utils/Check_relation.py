#------------- Explanation -------------#
"""
Author : Rick Vonk

This script can calculate the correlation bewteen two columns from csv_files

"""
#------------ Import -------------#
import numpy as np
import pandas as pd

#----------------- inputs ------------------#
csv_file = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/Combined_results/Swift_Combined_data.csv'

# ----------------- code ------------------ #

def df_maker(csv_file, column_wanted, column_id):
    df = pd.read_csv(csv_file)
    if isinstance(column_wanted, str):
        column_wanted = [column_wanted]
    
    columns_to_keep = [column_id] + list(column_wanted)
    df_modified = df[columns_to_keep]
    
    return df_modified

def df_averaging(df,column_id,  column_name, new_name):
    result = df.groupby(column_id)[column_name].mean().reset_index()
    result = result.rename(columns={column_name:new_name})
    
    return result

def correlation(df1, df2, sort_by):

    merged_df = pd.merge(df1, df2, on=sort_by, how='inner')
    col1 = [col for col in df1.columns if col != sort_by][0]
    col2 = [col for col in df2.columns if col != sort_by][0]

    correlation_matrix = merged_df[[col1, col2]].corr(method="spearman")

    relation_score = correlation_matrix.loc[col1, col2]

    print("--- Relationship Analysis ---")
    print(f"Correlation between {col1} and {col2}: {relation_score:.4f}")

    if relation_score > 0.7:
        print("Result: Strong positive relationship.")
    elif relation_score < -0.7:
        print("Result: Strong negative relationship.")
    elif abs(relation_score) < 0.3:
        print("Result: no relationship.")
    else:
        print("Result: Moderate relationship.")
    return


#----------------- Activation ------------------#
if __name__ == '__main__':

    df1_base = df_maker(csv_file, 'properly_docked_model_count', 'base_ID')
    #df1_avg = df_averaging(df_af_unf,'base_ID','AF3_confidence_score', 'AF3_score')

    df2_base = df_maker(csv_file,'avg_AF3_confidence','base_ID')
    #df2_avg = df_averaging(df_docking_unf,'base_ID','proper_docking','docking_factor')

    correlation(df1_base, df2_base, 'base_ID')


