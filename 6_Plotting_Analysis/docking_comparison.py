#------------- Explanation -------------#
"""
Author : Rick Vonk

This script is to check if there is a corrolation between true docking and proper orientation when using AF3

Input:
    - Total experiment folder

output:
    - csv-file with orientation data of all models
    - txt-file with analysis
    - boxplot with docking/Af3 correlation
    - a heatmap with dockign success rate

"""
#------------ Import -------------#
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#----------------- inputs ------------------#

Experiment_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

# ----------------- code ------------------ #

def Data_extraction(Exp_dir):
    df_list = []
    for case in glob.glob(os.path.join(Exp_dir,'*')):
        if 'tcr_alpha' not in case:
            continue
        data_file = os.path.join(case,'Intermediate_files','Confidence_docking_csv.csv')
        df_csv = pd.read_csv(glob.glob(data_file)[0])
        df_csv = df_csv.drop(columns=['cdr3_a','cdr3_b','tcr_a_seq','tcr_b_seq','Peptide'])
        df_list.append(df_csv)

    df_total = pd.concat(df_list, ignore_index=True)

    return df_total

def Data_modification(df_full):
    # expected format : experiment_id, AF3_confidence_score, docking_orientation, Reverse_true_docking
    # experiment_id : a[x]_b[x]_rs0_sample-[x]_model
    df_full['proper_docking'] = df_full['docking_orientation'] & df_full['Reverse_true_docking']
    df_full['base_ID'] = df_full['experiment_id'].str.extract(r'^(.*)_rs')
    df_full['a_group'] = df_full['experiment_id'].str.extract(r'(a\d+)')
    df_full['b_group'] = df_full['experiment_id'].str.extract(r'(b\d+)')

    df_modified = df_full.copy()
    # expected fromat returned : experiment_id,  AF3_confidence_score,  docking_orientation,  Reverse_true_docking,  proper_docking, base_ID, a_group, b_group
    return df_modified

def file_gen(output_dir):
    txt_file = os.path.join(output_dir, 'docking_analysis_file.txt')
    with open (txt_file,'w',newline='') as f:
        print('The following is the analysis of the docking orientations of the generated AF3 structures \n', file=f)
        print('Rick did a stupid and somehow removed all files for tcr_alpha_3', file=f)
    
    return txt_file

def docking_af3_corrolation_plot(df_modified, output_dir):
    fig, axes = plt.subplots(3, 1, figsize=(7, 12), sharey=True)
    
    df_modified.boxplot(column='AF3_confidence_score', by='docking_orientation', 
                        grid=False, ax=axes[0])
    axes[0].set_title('AF3 Confidence Score vs. Docking orientation, on top of MHC')
    axes[0].set_xlabel('Is TCR docked on top of MHC? (True = on-top / False = to the side)')
    axes[0].set_ylabel('AF3 Confidence Score')
    axes[0].text(-0.1, 1.05, 'A', transform=axes[0].transAxes, 
            fontsize=16, fontweight='bold', va='bottom', ha='right')

    df_modified.boxplot(column='AF3_confidence_score', by='Reverse_true_docking', 
                        grid=False, ax=axes[1])
    axes[1].set_title('AF3 Confidence Score vs. Reverse Docking orientation')
    axes[1].set_xlabel('Is properly docked or reversed? (True = proper / False = reversed)')
    axes[1].set_ylabel('AF3 Confidence Score')
    axes[1].text(-0.1, 1.05, 'B', transform=axes[1].transAxes, 
            fontsize=16, fontweight='bold', va='bottom', ha='right')

    df_modified.boxplot(column='AF3_confidence_score', by='proper_docking', 
                        grid=False, ax=axes[2])
    axes[2].set_title('AF3 Confidence Score vs. Proper Docking')
    axes[2].set_xlabel('Is Properly Docked? (True/False)')
    axes[2].set_ylabel('AF3 Confidence Score')
    axes[2].text(-0.1, 1.05, 'C', transform=axes[2].transAxes, 
            fontsize=16, fontweight='bold', va='bottom', ha='right')

    plt.suptitle(f'Plot of AF3-score and docking results, based on {len(df_modified.index)} models.',fontweight='bold') 
    
    plt.tight_layout()

    plot_out = os.path.join(output_dir, 'Correlation_plot_docking_AF3-score.png')
    plt.savefig(plot_out, dpi=300, bbox_inches='tight')
    plt.close()

    return

def Heatmap_success_rate(df_modified, output_dir):
    df = df_modified.copy()
    success_matrix = df.pivot_table(
    index="a_group", columns="b_group", values="proper_docking", aggfunc="mean"
    )

    sorted_b = sorted(success_matrix.columns, key=lambda x: int(x.replace("b", "")))
    sorted_a = sorted(success_matrix.index, key=lambda x: int(x.replace("a", "")))
    success_matrix = success_matrix.reindex(index=sorted_a, columns=sorted_b)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(success_matrix.values, cmap="YlGnBu", vmin=0, vmax=1, origin="upper")

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Success Rate", rotation=-90, va="bottom")
    cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    cbar.set_ticklabels(['(0/5)','(1/5)','(2/5)','(3/5)','(4/5)','(5/5)'])

    ax.set_xticks(np.arange(len(success_matrix.columns)))
    ax.set_yticks(np.arange(len(success_matrix.index)))
    ax.set_xticklabels(success_matrix.columns)
    ax.set_yticklabels(success_matrix.index)

    for i in range(len(success_matrix.index)):
        for j in range(len(success_matrix.columns)):
            val = success_matrix.values[i, j]
            if pd.isna(val):
                continue
            text_color = "white" if val > 1.0 else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontweight="bold",
            )

    ax.set_title("Proper Docking Success Rate per chain combination")
    ax.set_ylabel("A Chain")
    ax.set_xlabel("B Chain")
    plt.tight_layout()

    plot_out = os.path.join(output_dir, 'Heatmap_docking_success_rate.png')
    plt.savefig(plot_out, dpi=300, bbox_inches='tight')
    plt.close()
    return

def get_group_summary(dataframe, group_column):
    summary = dataframe.groupby(group_column).agg(
        mean_confidence=('AF3_confidence_score', 'mean'),
        proper_docking_rate=('proper_docking', 'mean'),
        orientation_rate=('docking_orientation', 'mean')
    ).reset_index()
    
    summary['%_proper_docking_rate'] = (summary['proper_docking_rate'] * 100).round(1)
    summary['%_orientation_rate'] = (summary['orientation_rate'] * 100).round(1)
    summary['mean_confidence'] = summary['mean_confidence'].round(2)
    
    summary = summary.sort_values(by='mean_confidence', ascending=False)

    final_columns = [group_column, 'mean_confidence', '%_proper_docking_rate', '%_orientation_rate']

    return summary[final_columns]

def DataFrame_counts(dataframe, column_name):
    Num_True = int(dataframe[column_name].sum())
    Num_False = int((~dataframe[column_name]).sum())
    total_count = len(dataframe.index)
    success_rate = round(Num_True/total_count,2)*100

    data = {
        'meaning': ['True', 'False', 'Success rate', 'Total'],
        column_name : [Num_True, Num_False, f'%{success_rate}', total_count]
    }
    
    return pd.DataFrame(data)

def Analysis(df_modified, output_dir, info_file):
    # top Af3 performing models 
    df_af3_sorted = df_modified.sort_values(by='AF3_confidence_score', ascending=False)
    top_per_id = df_af3_sorted.groupby('base_ID').head(1)

    top_per_id = top_per_id.drop(columns=[ 'docking_orientation',  'Reverse_true_docking',  'proper_docking', 'experiment_id', 'a_group', 'b_group'])
    Column_order = ['base_ID','AF3_confidence_score']

    top_per_id = top_per_id[Column_order].head(12)
    
    with open (info_file, 'a', newline='') as f:
        print('\n#----------------------------------------#', file=f)
        print('The top 12 performing models based on AF3-score by base_id', file=f)
        print(top_per_id.to_string(index=False), file=f)

    # look at AF3 score per chain
    a_group_summary = get_group_summary(df_modified, 'a_group')
    b_group_summary = get_group_summary(df_modified, 'b_group')

    with open (info_file, 'a', newline='') as f:
        print('\n#----------------------------------------#', file=f)
        print('A comparison of the AF3 scores per chain fro both alpha and beta\n', file=f)
        print("--- A-chain performance ---", file=f)
        print(a_group_summary.to_string(index=False), file=f)
        print("\n--- B-chain performance ---", file=f)
        print(b_group_summary.to_string(index=False), file=f)

    # Overall failure rates
    df_proper = DataFrame_counts(df_modified, 'proper_docking')
    df_orientation = DataFrame_counts(df_modified, 'docking_orientation')
    df_reverse = DataFrame_counts(df_modified, 'Reverse_true_docking')
    df_combined = df_proper.merge(df_orientation, on='meaning').merge(df_reverse, on='meaning')
    count = (~df_modified["docking_orientation"] & ~df_modified["Reverse_true_docking"]).sum()

    with open (info_file, 'a', newline='') as f:
        print('\n#----------------------------------------#', file=f)
        print('General info regarding true / false counts for the orientations\n', file=f)
        print(df_combined.to_string(index=False), file=f)
        print(f"The number of cases where both oriantation and reverse docking is False is {count}\n", file=f)

    return

def CSV_maker(df, result_dir, name):
    csv_file = os.path.join(result_dir, name)
    df.to_csv(csv_file)
    return


#----------------- Activation ------------------#
# Extra
Result_dir = os.path.join(Experiment_dir, 'Combined_results')
os.makedirs(Result_dir, exist_ok=True)

if __name__ == "__main__":
    Extracted_data = Data_extraction(Experiment_dir)
    modified_data = Data_modification(Extracted_data)

    txt_file = file_gen(Result_dir)

    docking_af3_corrolation_plot(modified_data, Result_dir)
    Heatmap_success_rate(modified_data, Result_dir)
    Analysis(modified_data, Result_dir, txt_file)

    CSV_maker(modified_data, Result_dir, 'combined_orientation_file.csv')