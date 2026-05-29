#------------- Explanation -------------#
"""
Author : Rick Vonk
This script is intended for after AF complex generation to extract information required for analysis.
This script will make a csv file containing structural and sequence info

Input :
    - An experiment name
    - Zipping statement

Output: 
    - if Zipping, than a tarred version of the 2_AF folder and zipping script used
    - an updates csv file with docking data
    - a csv-file with angles data
    - a tsv file for angle data
    - a txt-file with analytics

"""
#------------ Import -------------#
import os
import csv
import glob
import sys
import subprocess
import pandas as pd
from Bio import PDB
from pathlib import Path

#----------------- inputs ------------------#

experiment= 'tcr_alpha_11'

zipping = True

# ----------------- code ------------------ #
def tsv_file_maker(Working_dir):
    """
    This function makes a tsv file that is used in further functions, the tsv file contains chainn info for each model applied
    """
    files_dir = os.path.join(os.path.dirname(Working_dir), 'Intermediate_files')
    os.makedirs(files_dir,exist_ok=True)
    tsv_path = os.path.join(files_dir, 'tsv_for_orientation.tsv')
    pdb_loc = os.path.join(Working_dir, 'Process_2/TCR_post')
    with open(tsv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['pdb', 'Bchain', 'Achain', 'mhc_type', 'mhc_chain1', 'mhc_chain2', 'beta_organism', 'alpha_organism', 'mhc_chain1_organism', 'mhc_chain2_organism'])

        for entry in sorted(glob.glob(os.path.join(pdb_loc, '*.cif'))):
            name = os.path.splitext(os.path.basename(entry))[0]

            writer.writerow([
                name, 'B', 'A', 'MH1', 'M', 'N', 
                'homo sapiens', 'homo sapiens', 'homo sapiens', 'homo sapiens'
            ])
    return tsv_path

def cif2pdb(cif_filepath, pdb_filepath):
    """Convert a CIF file to PDB format using Biopython."""
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure('structure', cif_filepath)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_filepath)

def pdb_maker(Working_dir):
    """Iterates through cif folder and relocates and changed each file to a pdb(folder)"""
    cif_dir = os.path.join(Working_dir,'Process_2/TCR_post')
    cif_files = glob.glob(f'{cif_dir}/*.cif')
    output_dir = os.path.join(os.path.dirname(Working_dir),'input_pdb')
    os.makedirs(output_dir, exist_ok=True)
    for cif in cif_files:
        pdb_file = f"{output_dir}/{cif.split('/')[-1].replace('.cif', '.pdb')}"
        cif2pdb(cif, pdb_file)
    return output_dir

def crossing_angle_calc(pdb_loc, tsv_file, Working_dir, swiftTCR_utils,calc_incident_crossing_angle= '/home/rvonk1/3_Jolanda_data_Pipeline/Script_pipeline/2_csv-data/calc_incident_crossing_angle.py'):
    """
    This function uses calc_incident_crossing_angle.py to calculate the crossing angle of MHC
    """
    files_dir = os.path.join(os.path.dirname(Working_dir), 'Intermediate_files')

    subprocess.run([
        'python',
        calc_incident_crossing_angle,
        pdb_loc,
        tsv_file,
        swiftTCR_utils,
        files_dir
    ], check=True)

    CA_csv = os.path.join(files_dir,'angles_csv.csv')

    return CA_csv

def csv_addition(Working_dir, CA_csv):
    """Updates the csv file from the previous script with the crossing angle information"""
    files_dir = os.path.join(os.path.dirname(Working_dir), 'Intermediate_files')
    csv_file = os.path.join(files_dir, 'Confidence_docking_csv.csv')

    df_conf = pd.read_csv(csv_file)
    df_angles = pd.read_csv(CA_csv)

    df_angles['experiment_id'] = df_angles['modelname'].str.replace('.pdb', '', regex=False)
    angle_map = dict(zip(df_angles['experiment_id'], df_angles['crossing_angle']))
    df_conf['Reverse_true_docking'] = df_conf['experiment_id'].map(angle_map) > 0
    df_conf.to_csv(csv_file, index=False)

    return csv_file

def analytics_docking_crossing(df,variable_column, highlight):
    """Function used for analytics of the generated structures"""

    counts = df.groupby('ID')[variable_column].value_counts().unstack(fill_value=0)
    for col in [True, False]:
        if col not in counts:
            counts[col] = 0

    result_list = []
    for id_val, row in counts.iterrows():
        total = row[True] + row[False]
        perc_true = round((row[True] / total) * 100, 2) if total > 0 else 0
        
        output = f'{id_val} : {row[True]}, {row[False]}, {perc_true}%'
        if id_val in highlight:
            output += ' <--- "real"'
        result_list.append(f'{output}\n')
    return result_list

def run_analytics(csv_file,Working_dir, highlight):
    """
    Creates a .txt file with some analytics on the generated structures sorted by base_ID
    """

    files_dir = os.path.join(os.path.dirname(Working_dir), 'Intermediate_files')
    txt_file = os.path.join(files_dir,'csv_Analytics.txt')

    df = pd.read_csv(csv_file)
    id_col = df.columns[0]
    conf_col = df.columns[1]
    orient_col = 'docking_orientation'
    dock_col = 'Reverse_true_docking'

    parts = df[id_col].str.split('_')
    df['ID'] = parts.str[0] + "_" + parts.str[1]

    with open(txt_file,'w',newline='') as f:
        print('#----------- Confidence analytics ----------#', file=f)
        conf_stats = df.groupby('ID')[conf_col].agg(['mean', 'std']).round(2)

        print('ID : average, stdv', file=f)
        for id_val, row in conf_stats.iterrows():
            suffix = ' <--- "real"' if id_val in highlight else ''
            print(f'{id_val} : [{row["mean"]}, {row["std"]}]{suffix}', file=f)
    
        print('\n#----------- docking orientation analytics ----------#', file=f)
        print('ID : Total True, Total False, % True', file=f)
        print("".join(analytics_docking_crossing(df, orient_col, highlight)), file=f)
        
        print('\n#----------- reverse / normal docking analytics ----------#', file=f)
        print('ID : Total normal, Total reverse, % normal', file=f)
        print("".join(analytics_docking_crossing(df, dock_col, highlight)), file=f)

    return

def zipping_af3_files(template, Working_dir):
    """File to zip 2_AF files to reduce sproject space"""
    Working_dir = os.path.abspath(Working_dir)
    parent_dir = os.path.dirname(Working_dir)
    files_dir = os.path.join(parent_dir, 'Intermediate_files/zipping')
    os.makedirs(files_dir,exist_ok=True)

    with open(template, 'r') as f:
        contents = f.read()
    
    contents = contents.replace('$Result_location', files_dir)
    contents = contents.replace('$main_loc', parent_dir)
    if os.path.isfile(os.path.join(Working_dir,'2_AF.tar.gz')) :
        contents = contents.replace('$name', '2_AF_2')
    else:
        contents = contents.replace('$name', '2_AF')
    contents = contents.replace('$haddock_output_loc', Working_dir)

    script_path = os.path.join(files_dir, '2_AF_zipping')
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(contents)

    subprocess.run(['sbatch', script_path])

#----------------- Activation ------------------#
project_dir = os.path.join('/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments', experiment, '2_AF')
swiftTCR_utils = '/home/rvonk1/swifttcr/utils'

script_location = Path(__file__).resolve()
script_dir = script_location.parent
template_dir = os.path.join(script_dir,"required_for_scripts")
zip_template_loc = os.path.join(template_dir,'zipping_template.sh')


if __name__ == "__main__":
    tsv = tsv_file_maker(project_dir)
    pdb_location = pdb_maker(project_dir)
    CA_csv = crossing_angle_calc(
        pdb_location,
        tsv,
        project_dir,
        swiftTCR_utils,
        os.path.join(script_dir,'calc_incident_crossing_angle.py')
    )
    csv_file = csv_addition(project_dir, CA_csv)
    run_analytics(csv_file, project_dir,
        {'tcra2b3', 'tcra5b2', 'tcra1b1'}
    )

    if zipping:
        zipping_af3_files(
            zip_template_loc,
            project_dir
            )
