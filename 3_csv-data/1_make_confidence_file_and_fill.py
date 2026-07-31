#------------- Explanation -------------#
"""
Author : Rick Vonk
This script is intended for after AF complex generation to extract information required for analysis.
This script will make a csv file containing structural and sequence info

Input :
    - The Experiment directory, the overaching direction that is worked in.
    - Experiment name

Output :
    - Pipeline folders in the desired structure
    - csv file in 'Intermediate_files'with docking results

Usage :
    python /path/to/script/

"""
#------------ Import -------------#
import csv
import shutil
import os
import glob
import json
import re
import pandas as pd
from collections import defaultdict
from filter_docking_side import postfilter

#----------------- inputs ------------------#
Exp_dir = '/projects/0/prjs1135/report_Rick/4_Haddock_config_experimentation/experiments'

experiment= 'test_set_rs3'

# ----------------- code ------------------ #
def confidence_extraction(Working_dir):
    # Extracts AF3 confidence scores by mapping TCR_post files back to the AF3 inference output directory structure.
    confidence = {}
    search_path = os.path.join(Working_dir, 'Process_2/TCR_post', '*_rs*_sample-*_model.cif')
    
    for file_path in glob.glob(search_path):
        filename = os.path.basename(file_path)
        
        match = re.match(r'(.+)_rs(\d+)_sample-(\d+)(_model)', filename)            
        base_id, rs_val, sample_idx, suffix = match.groups()
        
        json_pattern = os.path.join(
            Working_dir, 
            'Process_1', 
            base_id, 
            'AF3_inference_output', 
            f'{base_id}_rs{rs_val}', 
            f'seed-*_sample-{sample_idx}', 
            'summary_confidences.json'
        )
        
        json_files = glob.glob(json_pattern)
        
        for jf in json_files:
            try:
                with open(jf, 'r') as f:
                    data = json.load(f)
                    experiment_key = f"{base_id}_rs{rs_val}_sample-{sample_idx}{suffix}"
                    confidence[experiment_key] = data.get('ranking_score')
            except (OSError, json.JSONDecodeError):
                continue

    return confidence
    
def csv_maker(Working_dir, confidence):
    # This function makes a csv file used later and fills it with the confidence score
    files_dir = os.path.join(os.path.dirname(Working_dir), 'Intermediate_files')
    os.makedirs(files_dir,exist_ok=True)
    csv_file = os.path.join(files_dir, 'Confidence_docking_csv.csv')
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['experiment_id', 'AF3_confidence_score','cdr3_a','cdr3_b','tcr_a_seq','tcr_b_seq','Peptide','docking_orientation','Reverse_true_docking'])
        for exp_id in sorted(confidence.keys()):
            writer.writerow([exp_id, confidence[exp_id],'','','','','',''])
    
    return csv_file

def aa3to1(res_name):
    # Convert three-letter amino acid code to one-letter code
    
    AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return AA3_TO_1.get(res_name.upper(), 'X')

def extract_residues(pdb_file):
    # Function used to extract residue info from a pdb file.
    seen = set()
    residues = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21].strip()
                res_seq = int(line[22:26])
                ins_code = line[26].strip()
                res_name = aa3to1(line[17:20].strip())

                key = (chain, res_seq, ins_code)

                if key not in seen:
                    seen.add(key)
                    residues.append((chain, res_seq, res_name))

    return residues

def cdr3_finder(Working_dir):
    # This function gives iterates over a pdb folder and extracts sequence information form each pdb file in the folder.
    seq_info = {}
    pdb_folder = os.path.join(Working_dir,'Process_2/Renumbered_TCR/renumbered')

    for file in os.listdir(pdb_folder):
        if not file.endswith(".pdb"): continue
        
        full_path = os.path.join(pdb_folder, file)
        info = extract_residues(full_path)
        ID = os.path.splitext(file)[0]

        tcr_a_seq_def = []
        tcr_b_seq_def = []
        cd3a = []
        cd3b = []

        for x in info :
            if x[0] == 'A':
                tcr_a_seq_def.append(x[2])
                # CDR3 in IMGT numbering lies between 105 and 118
                if x[1] in range(105, 118): 
                    cd3a.append(x[2])
            
            if x[0] == 'B':
                tcr_b_seq_def.append(x[2])
                if x[1] in range(105, 118):
                    cd3b.append(x[2])

        seq_info[ID] = [
            ''.join(cd3a),
            ''.join(cd3b),
            ''.join(tcr_a_seq_def),
            ''.join(tcr_b_seq_def),
            'RLSSCVPV'
            ]

    return seq_info

def csv_updater_cdr3(csv_file, data):
    # This function updates the csv made previously with the new information

    df = pd.read_csv(csv_file)
    df = df.astype(str).replace('nan', '')
    df = df.set_index('experiment_id')

    seq_df = pd.DataFrame.from_dict(
        data,
        orient='index',
        columns=['cdr3_a','cdr3_b','tcr_a_seq','tcr_b_seq','Peptide']
    )
    seq_df.index.name = 'experiment_id'
    df.update(seq_df)
    df = df.combine_first(seq_df)

    desired_order = [
        'AF3_confidence_score',
        'cdr3_a',
        'cdr3_b',
        'tcr_a_seq',
        'tcr_b_seq',
        'Peptide',
        'docking_orientation',
        'Reverse_true_docking'
    ]
    

    df = df.reindex(columns=desired_order)
    df.to_csv(csv_file)
    
    return

def docking_side_checker(csv_file, Working_dir):
    # This fucntion uses a fucntion from filter_docking_side.py to check if the CDR3 loop is within threshold distance of the peptide and returns True if yes, False if no.

    docking_side = {}

    with open(csv_file, 'r') as f:
        reading = csv.reader(f)
        header = next(reading)
        
        # Move the loop inside this block
        for case in reading:
            result = postfilter(
                f'{Working_dir}/Process_2/Renumbered_TCR/{case[0]}.pdb',
                case[6],
                case[2],
                case[3],
                min_num_residues=3,  
                threshold=15, 
                verbose=True,
                criteria='and', 
                tcr_dist_threshold=25
            )
            docking_side[case[0]] = result


    df = pd.read_csv(csv_file)
    df['docking_orientation'] = df['experiment_id'].map(docking_side)

    df.to_csv(csv_file, index=False)
    return

#----------------- Activation ------------------#

project_dir = os.path.join(Exp_dir, experiment, '2_AF')

if __name__ == "__main__":
    confidence_data = confidence_extraction(project_dir)
    csv_file = csv_maker(project_dir, confidence_data)
    info = cdr3_finder(project_dir)
    csv_updater_cdr3(csv_file, info)
    docking_side_checker(csv_file, project_dir)

