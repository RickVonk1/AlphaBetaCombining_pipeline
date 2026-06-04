#------------- Explanation -------------#
"""
Author : Rick Vonk

Extract chain information for SwiftTCR from existing structure data

Input:
    - The Experiment directory
    - An experiment name
    - Statement if zipping 2_AF

Output:
    - if zipping than a 2_AF.tar.gz file and removal of 2_AF folder
    - a directory with the SwiftTCR activation files
    - slurm jobs to activate and run SwiftTCR

"""
#------------ Import -------------#
import os
import glob
import subprocess
import pandas as pd
from pathlib import Path

#----------------- inputs ------------------#
exp_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'
SwiftTC_python_loc = '/home/rvonk1/swifttcr/scripts/swift_tcr.py'

Experiment_name = 'tcr_alpha_5'

zipping = False
# ----------------- code ------------------ #
def Swift_prepare(input_dir, swiftTCR_loc, template):
    input_structures = os.path.join(input_dir,'input_structures')

    for entry in os.listdir(input_dir):

        tcr_dir = os.path.join(input_dir, entry, 'input_structures','tcr')
        pmhc_dir = os.path.join(input_dir, entry, 'input_structures', 'pmhc')
        output = os.path.join(input_dir, entry, 'Swift_output')
        run_dir = os.path.join(input_dir, 'runs_dir')
        os.makedirs(output, exist_ok=True)
        os.makedirs(run_dir, exist_ok=True)

        with open (template, 'r') as f:
            content = f.read()
        
        content = content.replace('[case]', entry)
        content = content.replace('[output_dir]', output)
        content = content.replace('[tcr_loc]', tcr_dir)
        content = content.replace('[pmhc_loc]', pmhc_dir)

        with open (os.path.join(run_dir, f'Swift_run_{entry}.sh'), 'w') as f:
            f.write(content)
        
    return run_dir

def Swift_activation(run_dir):
    for entry in os.listdir(run_dir):
        if 'runs_dir' in entry:
            continue
            
        run_path = os.path.join(run_dir, entry) 
        command = ['sbatch', run_path]

        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)

        except subprocess.CalledProcessError as e:
            print(f"Failed to submit {entry}. Error: {e.stderr.strip()}")

def zipping_af3_files(template, Working_dir):
    """File to zip 2_AF files to reduce sproject space"""
    Working_dir = os.path.abspath(Working_dir)
    parent_dir = os.path.dirname(Working_dir)
    folder_name = os.path.basename(Working_dir)
    
    files_dir = os.path.join(parent_dir, 'Intermediate_files/zipping')
    os.makedirs(files_dir, exist_ok=True)

    with open(template, 'r') as f:
        contents = f.read()
    
    contents = contents.replace('$Result_location', files_dir)
    contents = contents.replace('$main_loc', parent_dir)
    contents = contents.replace('$folder_to_zip', folder_name)
    
    if os.path.isfile(os.path.join(Working_dir, '2_AF.tar.gz')):
        contents = contents.replace('$name', '2_AF_2')
    else:
        contents = contents.replace('$name', '2_AF')
        
    contents = contents.replace('$haddock_output_loc', Working_dir)

    script_path = os.path.join(files_dir, '2_AF_zipping')
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(contents)

    subprocess.run(['sbatch', script_path])
    return
#----------------- Activation ------------------#

swift_dir = os.path.join(exp_dir, Experiment_name, 'SwiftTCR')

script_location = Path(__file__).resolve()
script_dir = script_location.parent
template = os.path.join(script_dir, 'Template','template_run_swifttcr.sh')
zip_template_loc = os.path.join(script_dir.parent,'4_haddock',"required_for_scripts",'zipping_template.sh') 


if __name__ == "__main__":

    runs_dir = Swift_prepare(swift_dir, SwiftTC_python_loc, template)
    Swift_activation(runs_dir)
    
    if zipping:
        zipping_af3_files(
            zip_template_loc,
            os.path.join(exp_dir, Experiment_name, '2_AF')
            )
