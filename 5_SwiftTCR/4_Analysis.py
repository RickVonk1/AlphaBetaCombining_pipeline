#------------- Explanation -------------#
"""
Author : Rick Vonk

Extract chain information for SwiftTCR from existing structure data

Input:
    -

Output:
    -

"""
#------------ Import -------------#
import os
import glob
import re
import pandas as pd
import math
import subprocess
from pathlib import Path


#----------------- inputs ------------------#
exp_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

Experiment_name = 'tcr_alpha_1'

zipping = False

# ----------------- code ------------------ #
def compute_entropy(sizes):
    total_elements = sum(sizes)
    if total_elements == 0:
        return 0.0, 0
    
    entropy = 0.0
    for size in sizes:
        if size > 0:
            p = size / total_elements
            entropy -= p * math.log2(p)
            
    return entropy, total_elements

def calculate_shannon_entropy(experiment_dir):
    Swift_output = os.path.join(experiment_dir, 'SwiftTCR')
    Cluster_result = []
    for case in glob.glob(os.path.join(Swift_output,'*')):
        if 'runs_dir' in case:
            continue

        base_ID = os.path.basename(case)
        combined_file = os.path.join(case, 'Swift_combined', base_ID[:4], 'clustering.txt')
        sizes_just_neighbors = []
        pattern = re.compile(r'with (\d+) neighbors')
    
        with open(combined_file, 'r') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    neighbors = int(match.group(1))
                    sizes_just_neighbors.append(neighbors + 1)

        entropy_neighbors, total_n = compute_entropy(sizes_just_neighbors)
        Cluster_result.append({'base_ID': base_ID, 'Shannon_entropy_score': entropy_neighbors})

    df = pd.DataFrame(Cluster_result)
    
    return df

def Swift_extraction(experiment_dir):
    Swift_output = os.path.join(experiment_dir, 'SwiftTCR')
    Cluster_result = []
    for case in glob.glob(os.path.join(Swift_output,'*')):
        if 'runs_dir' in case:
            continue

        base_ID = os.path.basename(case)
        combined_file = os.path.join(case, 'Swift_combined', base_ID[:4], 'clustering.txt')
        with open (combined_file, 'r') as f:
            content = f.readline()
            match = re.search(r'with\s+(\d+)\s+neighbors', content)
            top_cluster = int(match.group(1))
            Cluster_result.append({'base_ID': base_ID, 'Top_Cluster_size': top_cluster})
    df = pd.DataFrame(Cluster_result)

    return df

def confidence_extraction(experiment_dir):
    confidence = os.path.join(experiment_dir, 'Intermediate_files', 'Confidence_docking_csv.csv')
    df = pd.read_csv(confidence)
    
    df_useful = df[['experiment_id', 'docking_orientation', 'Reverse_true_docking', 'AF3_confidence_score']].copy()
    
    df_useful['base_ID'] = df_useful['experiment_id'].str.split('_rs').str[0]

    proper_docked = (df_useful['docking_orientation'] == True) & (df_useful['Reverse_true_docking'] == True)
    df_filtered = df_useful[proper_docked]
    
    proper_base_ID = df_filtered.groupby('base_ID').agg(
        properly_docked_model_count=('experiment_id', 'size'),
        avg_AF3_confidence=('AF3_confidence_score', 'mean')
    ).reset_index()

    return proper_base_ID

def csv_maker(experiment_dir):
    df_combined_unf = pd.merge(Swift_extraction(experiment_dir), confidence_extraction(experiment_dir), on='base_ID')
    df_combined = pd.merge(df_combined_unf, calculate_shannon_entropy(experiment_dir), on='base_ID').round(3)
    csv_loc = os.path.join(experiment_dir,'Intermediate_files','Swift_Result.csv')

    Swift_csv = df_combined.to_csv(csv_loc,index=False)
    return 

def zipping_af3_files(template, Working_dir):
    """File to zip SwiftTCR files to reduce sproject space"""
    Working_dir = os.path.abspath(Working_dir)
    parent_dir = os.path.dirname(Working_dir)
    folder_name = os.path.basename(Working_dir)
    
    files_dir = os.path.join(parent_dir, 'Intermediate_files',folder_name,)
    os.makedirs(files_dir, exist_ok=True)

    with open(template, 'r') as f:
        contents = f.read()
    
    contents = contents.replace('$Result_location', files_dir)
    contents = contents.replace('$main_loc', parent_dir)
    contents = contents.replace('$folder_to_zip', folder_name)
    
    if os.path.isfile(os.path.join(Working_dir, f'{folder_name}.tar.gz')):
        contents = contents.replace('$name', f'{folder_name}_2')
    else:
        contents = contents.replace('$name', f'{folder_name}')
        
    contents = contents.replace('$haddock_output_loc', Working_dir)

    script_path = os.path.join(files_dir, f'{folder_name}_AF_zipping')
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(contents)

    subprocess.run(['sbatch', script_path])
    return
#----------------- Activation ------------------#
script_location = Path(__file__).resolve()
script_dir = script_location.parent
zip_template_loc = os.path.join(script_dir.parent,'4_haddock',"required_for_scripts",'zipping_template.sh') 


experiment_dir = os.path.join(exp_dir, Experiment_name)

if __name__ == "__main__":
    
    csv_maker(experiment_dir)
    
    if zipping:
        zipping_af3_files(
            zip_template_loc,
            os.path.join(exp_dir, Experiment_name, 'SwiftTCR')
            )