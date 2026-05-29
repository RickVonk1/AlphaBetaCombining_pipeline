import sys
import os
import glob
import warnings
import subprocess
from Bio.PDB import MMCIFParser, PDBIO
from argparse import ArgumentParser

"""
Script converts mmcif to pdb and renumbers them using ANARCI's immunopdb script.

ANARCI github: https://github.com/oxpig/ANARCI
"""

def main(mmcif_path,outputdir, immunopdb):
    """
    This script converts a mmcif to a pdb file and renumbers the TCR structure
    to IMGT numbering appying ANARCI.
    """
    output_pdb = mmcif_to_pdb(mmcif_path,outputdir)
    run_anarci(output_pdb, outputdir, immunopdb)

def run_command(command):
    """
    Helper function to run a command with subprocess.run and handle errors.
    
    Args:
        command (str): The shell command to be executed.
        
    Returns:
        str: The standard output of the command if successful.
    
    Raises:
        subprocess.CalledProcessError: If the command fails.
    """

    #Using os.popen provides a clearer error logging
    command = (' ').join(command)
    result = os.popen(command).read()
    return result

    # try:
        #result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        #return result.stdout
    # except subprocess.CalledProcessError as e:
    #     print(f"Error occurred while running command: {command}")
    #     print(f"Error message: {e.stderr}")
    #     raise


def mmcif_to_pdb(mmcif_path, output_dir):
    """
    This function converts mmcif files to pdb. 

    input:  mmcif_path = path to the mmcif file that needs to be converted to pdb.
            output_dir = the directory to place the pdb file.
    
    output: A pdb will be saved in the directory containing the same
            basename as the mmcif and the path will be returned. 
    
    """
    basename = os.path.basename(mmcif_path).split(".")[0]
    pdb_file = f"{basename}.pdb"

    parser = MMCIFParser()
    structure = parser.get_structure("model", mmcif_path)

    output_pdb = os.path.join(output_dir,pdb_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

    return output_pdb

def fix_TER_newline(filepath):
    '''Fixes the formatting of a PDB file so that TER and END records are on new lines'''
    with open(filepath, 'r') as infile:
        filedata = infile.read()
        filedata = filedata.replace('              ATOM', '            \nATOM')
        filedata = filedata.replace('              END', '            \nEND')
    with open(filepath, 'w') as outfile:
        outfile.write(filedata)

def run_anarci(pdb_path,output_dir, immunopdb):
    #shift TCR numbering due to overwriting residue numbers
    #to extract only TCRa and TCRb (A and B) chains, you can use:
    #  pdb_selchain -A {pdb_path} > tcra.pdb 
    #  pdb_selchain -B {pdb_path} > tcrb.pdb
    #  pdb_merge tcra.pdb tcrb.pdb > tcr_only.pdb
    basename = os.path.basename(pdb_path).replace(".pdb","")
    output_path = os.path.join(output_dir, 'renumbered', basename +'.pdb')

    temp_file = f"{output_dir}/{basename}_tcr_shifted.pdb"
    command_sel_tcr_reres = f"pdb_reres -500 {pdb_path} > {temp_file}"
    #run_command(command_sel_tcr_reres)
    #command_sel_tcr_reres = (' ').join(command_sel_tcr_reres)
    os.popen(command_sel_tcr_reres).read()

    if not os.path.exists(temp_file):
        raise FileNotFoundError(f"Error: The expected output file {temp_file} was not created.")

    fix_TER_newline(temp_file)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", module='Bio.PDB')
        command = [
            "python", immunopdb,
            "-i", temp_file,
            "-o", output_path,
            "-s", "imgt",
            "--receptor", "tr"
        ]
        command = (' ').join(command)
        result = os.popen(command).read()
        # subprocess.run([
        #     "python", immunopdb,
        #     "-i", f"{output_dir}/{basename}_tcr_shifted.pdb",
        #     "-o", output_path,
        #     "-s", "imgt",
        #     "--receptor", "tr"
        # ], check=True)
    if not os.path.exists(output_path):
        print(f"Error: The expected output file {output_path} was not created after running ANARCI.")
        print(f"Commad used: {command}")
        raise FileNotFoundError(f"Error: ANARCI did not produce the expected output file {output_path}.")

    try:
        os.remove(temp_file)
    except FileNotFoundError:
        print(f"Warning: Temporary file {temp_file} not found for deletion.")

if __name__ == "__main__":

    parser = ArgumentParser(description="Renumbers TCR with ANARCI's immunopdb and outputs them as pdbs")
    
    parser.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory containing RMSD .txt files from ProFit",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory where clustering outputs will be written",
    )
    parser.add_argument(
        "--immunopdb-path",
        type=str,
        help="Path to the ImmunoPDB script of ANARCI",
        default="/projects/0/prjs1135/software/ANARCI/Example_scripts_and_sequences/ImmunoPDB.py"
    )

    args = parser.parse_args()
    immunopdb = args.immunopdb_path
    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(output_dir + "/renumbered"):
        os.makedirs(output_dir + "/renumbered")

    #loop over input dir 
    for path in glob.glob(f"{input_dir}/*"):
        id_output = os.path.basename(path).split(".")[0]
        if  path.endswith(".cif"):
            main(path, output_dir, immunopdb)
