#------------- Explanation -------------#
"""
Author : Rick Vonk

This is the first script in the haddock3 pipeline, this script is able to make an experiment directory with required subdirectories and able to put relavant pdb entries in the input folder

Input :
    - The Experiment directory
    - An experiment name

Output:
    - Remapped input structures
    - Folders for future use in the pipeline

"""
#------------ Import -------------#
import csv
import shutil
import os
import pandas as pd
import glob
from Bio import PDB

#----------------- inputs ------------------#
working_dir = '/projects/0/prjs1135/report_Rick/4_Haddock_config_experimentation/experiments'

experiment = 'test_set_redo'

# ----------------- code ------------------ #
def Folder_maker(experiment_name, output_dir):   
    exp_path = os.path.join(output_dir, experiment_name)
    pdb_loc = os.path.join(exp_path, "input_pdb")
    result_loc = os.path.join(exp_path, "Result")
    haddock_out = os.path.join(exp_path, "Haddock_output")
    config_dir = os.path.join(haddock_out, "config_files")
    os.makedirs(exp_path,exist_ok=True)
    os.makedirs(pdb_loc, exist_ok=True)
    os.makedirs(result_loc,exist_ok=True)
    os.makedirs(config_dir, exist_ok=True)
    return pdb_loc

def pdb_input_remapping(source_folder, orientation_csv):
    df = pd.read_csv(orientation_csv)
    valid_ids = df[(df['docking_orientation'] == True) & (df['Reverse_true_docking'] == True)]['experiment_id'].tolist()
    count = 0
    
    parser = PDB.PDBParser(QUIET=True)
    io = PDB.PDBIO()
    
    for file_name in os.listdir(source_folder):
        if file_name.endswith(".pdb"):
            exp_id = os.path.splitext(file_name)[0]            
            if exp_id in valid_ids:
                file_path = os.path.join(source_folder, file_name)
                structure = parser.get_structure('protein', file_path)
                mapping = {'B': 'A', 'N': 'M'}

                for model in structure:
                    for old_id, new_id in mapping.items():
                        
                        if old_id in model and new_id in model:
                            target_chain = model[new_id]
                            source_chain = model[old_id]
                            
                            res_nums = [res.id[1] for res in target_chain.get_residues()]
                            offset = (max(res_nums) if res_nums else 0) + 1000

                            for residue in list(source_chain):
                                source_chain.detach_child(residue.id)
                                old_res_id = residue.id
                                new_res_id = (old_res_id[0], old_res_id[1] + offset, old_res_id[2])
                                residue.id = new_res_id                                
                                target_chain.add(residue)
                            
                            model.detach_child(old_id)
                        
                        elif old_id in model:
                            source_chain = model[old_id]
                            
                            for residue in list(source_chain):
                                source_chain.detach_child(residue.id)
                                old_res_id = residue.id
                                new_res_id = (old_res_id[0], old_res_id[1] + 1000, old_res_id[2])
                                residue.id = new_res_id
                                source_chain.add(residue)
                                

                            source_chain.id = new_id

                io.set_structure(structure)
                io.save(file_path)
                count += 1
    return

#----------------- Activation ------------------#
Orientation = os.path.join(working_dir,experiment,'Intermediate_files', 'Confidence_docking_csv.csv')

if __name__ == "__main__":
    pdb_loc = Folder_maker(experiment, working_dir)
    pdb_input_remapping(pdb_loc, Orientation)




