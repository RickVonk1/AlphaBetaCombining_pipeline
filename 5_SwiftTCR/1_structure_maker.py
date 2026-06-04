#------------- Explanation -------------#
"""
Author : Rick Vonk

Extract chain information for SwiftTCR from existing structure data

Input:
    - The Experiment directory
    - An experiment name

Output:
    - directory with input structures

"""
#------------ Import -------------#
import os
import glob
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select

#----------------- inputs ------------------#
exp_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'

Experiment_name = 'tcr_alpha_1'

# ----------------- code ------------------ #
class ChainSelect(Select):
    def __init__(self, chains_to_keep):
        self.chains_to_keep = chains_to_keep

    def accept_chain(self, chain):
        return 1 if chain.get_id() in self.chains_to_keep else 0

def split_pdb(input_file, tcr_chains, pmhc_chains, sample_name, base_name, output):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tcr_pmhc", input_file)
    
    io = PDBIO()
    io.set_structure(structure)
    
    tcr_dir = os.path.join(output, base_name,'input_structures','tcr')
    os.makedirs(tcr_dir, exist_ok=True)
    tcr_pdb = os.path.join(tcr_dir,f'{base_name}_tcr_{sample_name}.pdb')
    io.save(tcr_pdb, select=ChainSelect(tcr_chains))
    
    pmhc_dir = os.path.join(output, base_name,'input_structures','pmhc')
    os.makedirs(pmhc_dir, exist_ok=True)
    pmhc_pdb = os.path.join(pmhc_dir,f'{base_name}_pmhc_{sample_name}.pdb')
    io.save(pmhc_pdb, select=ChainSelect(pmhc_chains))

    return

def input_sctruc(experiment_dir, experiment_name, output_dir, orientation_csv):
    pdb_loc = os.path.join(experiment_dir, experiment_name,'2_AF','Process_2','Renumbered_TCR', 'renumbered')

    df = pd.read_csv(orientation_csv)
    valid_ids = df[(df['docking_orientation'] == True) & (df['Reverse_true_docking'] == True)]['experiment_id'].tolist()
    
    for pdb in glob.glob(os.path.join(pdb_loc,'*')):
        base_name = os.path.basename(pdb).split('.pdb')[0]
        if base_name not in valid_ids:
            continue

        parts = base_name.split('_')
        base_ID = f'{parts[0]}_{parts[1]}'
        sample_ID = parts[3]

        split_pdb(
            input_file=pdb, 
            tcr_chains=['A','B'],
            pmhc_chains=['M','N','P'],
            sample_name= sample_ID,
            base_name= base_ID,
            output= output_dir
        )

    return

#----------------- Activation ------------------#

swift_dir = os.path.join(exp_dir, Experiment_name, 'SwiftTCR')
orientation = os.path.join(exp_dir, Experiment_name, 'Intermediate_files','Confidence_docking_csv.csv')

if __name__ == "__main__":
    input_sctruc(exp_dir, Experiment_name, swift_dir, orientation)